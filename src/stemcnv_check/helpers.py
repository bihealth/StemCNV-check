# -*- coding: utf-8 -*-
"""Helper functions for pipeline"""
import importlib.resources
import pandas as pd
import os
import re
import ruamel.yaml as ruamel_yaml
from pathlib import Path
from pydantic.v1.utils import deep_update
from . import STEM_CNV_CHECK, MEHARI_DB_VERSION, ENSEMBL_RELEASE
from .exceptions import SampleConstraintError, ConfigValueError, CacheUnavailableError
from collections import OrderedDict
from loguru import logger as logging


def read_sample_table(filename, name_remove_regex=None, return_type='dataframe'):
    """Read sample table from file, return either a list of the required col entires,
    a list with Ordereddict of all rows (colnames as keys) or a pandas dataframe"""

    if return_type not in ('list', 'list_withopt', 'dataframe'):
        raise ValueError('Invalid return_type: ' + return_type)

    if str(filename).endswith('.xlsx'):
        def reader_func(f): return pd.read_excel(f, comment='#', dtype=str).dropna(how='all')
    elif str(filename).endswith('.tsv') or str(filename).endswith('.txt'):
        def reader_func(f): return pd.read_csv(f, sep='\t', comment='#')
    elif str(filename).endswith('.csv'):
        def reader_func(f): return pd.read_csv(f, comment='#')
    else:
        raise ValueError('Unknown file format for sample table: ' + filename)

    req_cols = [
        'Sample_ID', 'Chip_Name', 'Chip_Pos', 'Array_Name', 'Sex',
        'Reference_Sample', 'Regions_of_Interest', 'Sample_Group'
    ]

    if name_remove_regex:
        logging.debug(f'Removing regex from column names: "{name_remove_regex}"')
        def rename_func(name): return re.sub(name_remove_regex, '', name)
    else:
        def rename_func(name): return name

    sample_tb = reader_func(filename).rename(columns=rename_func, errors='raise').fillna('').astype(str)
    # logging.debug('Final column names:' + ', '.join(sample_tb.columns))
    # FIXME: raise error on duplicated col names?

    if not all(col in sample_tb.columns for col in req_cols):
        logging.debug('Sample table columns: ' + ', '.join(sample_tb.columns))
        logging.debug('Required columns: ' + ', '.join(req_cols))

        missing = [c for c in req_cols if c not in sample_tb.columns]
        raise SampleConstraintError(
            'Not all required sample_table columns found. Missing columns: ' + ', '.join(missing)
        )
    # reorder columns to required order
    cols = req_cols + [col for col in sample_tb.columns if col not in req_cols]
    sample_tb = sample_tb[cols]

    if return_type == 'list':
        sample_list = []
        for _, row in sample_tb.iterrows():
            sample_list.append([row[col] for col in req_cols])
        return sample_list
    elif return_type == 'list_withopt':
        sample_list = []
        for _, row in sample_tb.iterrows():
            sample_list.append(OrderedDict((col, row[col]) for col in sample_tb.columns))
        return sample_list
    elif return_type == 'dataframe':
        return sample_tb.set_index('Sample_ID', drop=False)


def make_apptainer_args(config, cache_path, tmpdir=None, not_existing_ok=False, extra_bind_args=None):
    """Collect all outside filepaths that need to be bound inside container"""

    bind_points = [
        (config['data_path'], '/outside/data'),
        (config['raw_data_folder'], '/outside/rawdata'),
        (config['log_path'], '/outside/logs'),
        (importlib.resources.files(STEM_CNV_CHECK), '/outside/snakedir'),
    ]
    # When apptainer is used the global files should always be present
    used_genomes = set(array['genome_version'] for name, array in config['array_definition'].items())
    for global_file in ('fasta', 'gtf', 'genome_info', 'mehari_txdb', 'dosage_scores'):
        for genome_version in used_genomes:
            static_file = get_global_file(global_file, genome_version, config['global_settings'], cache_path)
            bind_points.append((static_file, '/outside/static/' + os.path.basename(static_file)))

    if tmpdir is not None:
        bind_points.append((tmpdir, '/outside/tmp'))

    # Apptainer is needed to create some of the array definition files, in these cases "not_existing_ok" should be set to True
    # However, the directory in which the files should be created still need to be linked
    array_base_dirs = {}
    for array in config['array_definition'].keys():
        for name, file in config['array_definition'][array].items():
            if name == 'genome_version':
                continue
            # CSV file is optional
            elif name == 'csv_manifest_file' and not file:
                continue
            # Replace default cache path with actual cache path
            if (name.startswith('penncnv') or name.startswith('array')) and file == '__cache-default__':
                file = get_array_file(name, array, config, cache_path)
            # Can only directly mount existing files
            if not os.path.isfile(file):
                # Only allow missing files for penncnv and array files
                if not_existing_ok and (name.startswith('penncnv') or name.startswith('array')):
                    if file == '__cache-default__':
                        file = get_array_file(name, array, config, cache_path)
                    # Try to bind the base directory where files should be written, if it exists
                    # Note: If multiple files from one array have different base directories, this will fail
                    file_basedir = os.path.dirname(file)
                    # If the filename not a path (i.e. string/in the local directory) the basedir needs to be set to pwd
                    if not file_basedir:
                        file_basedir = os.getcwd()
                    if not os.path.isdir(file_basedir):
                        logging.info(f"Creating base directory for array static-data '{array}': {file_basedir}")
                        os.makedirs(file_basedir)
                    dir_mount = (
                        file_basedir.format(cache=cache_path, array_name=array),
                        f'/outside/writearray/{array}'
                    )
                    if array in array_base_dirs:
                        if array_base_dirs[array] != file_basedir:
                            logging.error(
                                f"Multiple base directories found for static data files of array '{array}': "
                                f"{array_base_dirs[array]} and {file_basedir}. Static-data workflow "
                                f"will likely to fail writing to the intended file paths."
                            )
                    else:
                        array_base_dirs[array] = file_basedir
                        bind_points.append(dir_mount)
                    continue
                else:
                    raise FileNotFoundError(f"Array definition file '{file}' does not exist  ({name} for {array}).")
            bind_points.append((file, '/outside/{array}/{fname}'.format(array=array, fname=os.path.basename(file))))

    # Sort bind_points to make testing easier (loops over config_dict add them in a non-deterministic order)
    # Also convert them to Path objects and resolve and links (which might not work once inside the container)
    bind_point_str_list = []
    for host, cont in sorted(bind_points, key=lambda x: x[1]):
        host = Path(host).resolve()
        # Double check that all to-be-mounted files exist
        if not host.exists():
            raise FileNotFoundError(
                f"File or directory '{host}' does not exist, but was supposed to be mounted for use with apptainer."
            )
        bind_point_str_list.append(f"'{host}':'{cont}'")
    bind_point_str = "-B " + ','.join(bind_point_str_list)

    # allow additional arguments to be passed
    if extra_bind_args:
        if isinstance(extra_bind_args, str):
            extra_bind_args = [extra_bind_args]
        bind_point_str = bind_point_str + ',' + ','.join(list(extra_bind_args))

    logging.debug("Binding points for apptainer: " + str(bind_point_str))

    return bind_point_str


def collect_SNP_cluster_ids(sample_id, clustering_config, sample_data_df):

    def check_collected_ids(add_ids, max_n):
        # remove the original sample_id
        if sample_id in add_ids:
            add_ids.remove(sample_id)
        # remove duplicates, but keep order
        add_ids = list(OrderedDict.fromkeys(add_ids))

        # If any collected ids don't exist, fail
        if not set(add_ids).issubset(sample_data_df.index):
            not_found = set(add_ids) - set(sample_data_df.index)
            raise SampleConstraintError(
                "Some of the extracted sample_id's SNP clustering do not exist in the sample_table: " + ', '.join(not_found)
            )
        # Check that all samples belong to the same array, remove all that don't match
        sample_array = sample_data_df['Array_Name'].loc[sample_id]
        wrong_array = set()
        for id in add_ids:
            if sample_data_df['Array_Name'].loc[id] != sample_array:
                wrong_array.add(id)
        if wrong_array:
            logging.warning(
                f"Samples {', '.join(wrong_array)} do not belong to the same array as {sample_id} ({sample_array}). "
                f"They will be excluded from SNP clustering."
            )
            add_ids = [id for id in add_ids if id not in wrong_array]

        # Check that the number of samples is not too large, only count _new_ additions
        n_postmerge = len(ids | set(add_ids))
        skip_ids = []
        if max_n and len(ids) >= max_n:
            skip_ids = add_ids[:]
            add_ids = []
        elif max_n and n_postmerge > max_n:
            max_add = max_n - len(ids)
            add_ids, skip_ids = add_ids[:max_add], add_ids[max_add:]

        if skip_ids:
            logging.warning(
                f"Too many samples for SNP clustering of {sample_id} ({n_postmerge}), only the first {max_n} will be used.",
                once=True
            )
            logging.warning(f"Skipping: {', '.join(skip_ids)}")

        return add_ids

    ids = set()
    # Add reference sample if it exists
    if sample_data_df['Reference_Sample'].loc[sample_id]:
        ids.add(sample_data_df['Reference_Sample'].loc[sample_id])

    # sample_ids as is
    checked_ids = check_collected_ids(clustering_config['sample_ids'], clustering_config['max_number_samples'])
    ids.update(checked_ids)

    # 'id_columns' entry: take all sample_ids from '[column]'
    for col in clustering_config['id_columns']:
        if col not in sample_data_df.columns:
            raise ConfigValueError('Config for SNP clustering refers to non-existing column: ' + col)
        checked_ids = check_collected_ids(
            [id for id in sample_data_df[col].loc[sample_id].split(',') if id],
            clustering_config['max_number_samples']
        )
        ids.update(checked_ids)

    # 'match_columns' entries: take all sample_ids with the same value in '[column]'
    for col in clustering_config['match_columns']:
        if col not in sample_data_df.columns:
            raise ConfigValueError('Config for SNP clustering refers to non-existing column: ' + col)
        match_val = sample_data_df[col].loc[sample_id]
        # Only actually use truthy values (not empty strings or NAs)
        if match_val:
            checked_ids = check_collected_ids(
                sample_data_df.set_index(col)['Sample_ID'].loc[[match_val]].to_list(),
                clustering_config['max_number_samples']
            )
            ids.update(checked_ids)

    return ids


def config_extract(entry_kws, config, def_config):
    key_tree = list(entry_kws)
    subconfig = config
    subconfig_def = def_config
    used_entries = []
    while key_tree:
        entry = key_tree.pop(0)
        used_entries.append(entry)
        if entry in subconfig:
            subconfig = subconfig[entry]
            subconfig_def = subconfig_def[entry]
        elif entry in subconfig_def:
            subconfig = subconfig_def[entry]
            subconfig_def = subconfig_def[entry]
            logging.debug('Using config default values for: ' + ' : '.join(used_entries))
        else:
            logging.warning('"' + ' : '.join(used_entries) + '" is not a valid config entry or has been deprecated')
            # warnings.warn(' : '.join(used_entries) + " is not a valid config entry or has been deprecated",
            #               ConfigKeyWarning)
            return None

    return subconfig


def load_config(args, inbuilt_defaults=True):
    yaml = ruamel_yaml.YAML(typ='safe')
    with open(args.config) as f:
        config = yaml.load(f)

    # Use the global array definition block, if it exists: update with user-defined values
    if not args.no_cache:
        # check whether a global array definition file exists in the cache and load it
        # if it doesn't exist or doesn't contain the array definition block, return the config as is
        cache_array_defs = get_cache_array_definition(get_cache_dir(args, config))
        if cache_array_defs:
            with open(cache_array_defs, 'r') as f:
                global_array_config = yaml.load(f)

            if 'array_definition' not in global_array_config:
                pass
            # Merge the configs, with the user-defined values taking precedence
            elif 'array_definition' not in config or not config['array_definition']:
                config['array_definition'] = global_array_config['array_definition']
            else:
                config['array_definition'] = deep_update(global_array_config['array_definition'], config['array_definition'])

    # Update inbuilt default config with all previously defined values
    if inbuilt_defaults:
        config_def_file = importlib.resources.files(STEM_CNV_CHECK).joinpath('control_files', 'default_config.yaml')
        with config_def_file.open() as f:
            default_config = yaml.load(f)
        config = deep_update(default_config, config)

    return config

def get_default_cache_path():
    yaml = ruamel_yaml.YAML(typ='safe')
    with importlib.resources.files(STEM_CNV_CHECK).joinpath('control_files').joinpath('default_config.yaml').open() as f:
        default_config = yaml.load(f)
    return Path(default_config['global_settings']['cache_dir']).expanduser()


def get_cache_dir(args, config):
    """Check if given path can be used for caching, get default (in install dir or home) otherwise, create the path"""
    def check_path_usable(path):
        try:
            if path.is_dir() and os.access(path, os.W_OK):
                logging.debug(f"Using existing cache directory: {path}")
            else:
                path.mkdir(parents=True)
                logging.debug(f"Created cache directory: {path}")
            return path
        except PermissionError:
            logging.debug(f"Failed to create or access cache in: {path}")
        except FileExistsError:
            logging.debug(f"Cache path exists but is not accessible: {path}")
        return None
    if args.no_cache:
        logging.info("No cache directory will be used, conda and docker images will be stored in snakemake project directory")
        return None
    elif args.cache_path:
        # Check if user supplied path is usable
        cache_path = check_path_usable(Path(args.cache_path))
        if cache_path:
            return cache_path
        else:
            logging.debug(f"Could not use cache directory '{args.cache_path}'")
    else:
        auto_paths = []
        if 'global_settings' in config and 'cache_dir' in config['global_settings']:
            auto_paths.append(('config', Path(config['global_settings']['cache_dir']).expanduser()))
        auto_paths.append(('home', get_default_cache_path()))

        for name, path in auto_paths:
            cache_path = check_path_usable(path)
            if cache_path:
                return cache_path
            else:
                logging.warning(f"Could not use cache directory '{path}'")

    # Nothing worked, return None & don't use specific cache
    logging.warning("Cache directory is not usable! Trying to run without cache. Conda and docker images will be stored in snakemake project directory")
    return None

def get_cache_array_definition(cache_path, config_file='global_array_definitions.yaml'):
    """Get config file with global array definitions from cache path if it exists"""
    if cache_path:
        global_config = os.path.join(cache_path, config_file)
        if os.path.isfile(global_config):
            return global_config
    return ''


def get_global_file(type, genome_version, global_settings, cache_path, fill_wildcards=True):
    """Get path to files defined in the global settings. All of these have an internally defined default.
    Supported types are: fasta, gtf, genome_info, mehari_txdb"""

    genome = 'hg38' if genome_version in ("hg38", "GRCh38") else 'hg19'
    wildcards = {'genome': genome_version}

    if type == 'fasta':
        config_entry = global_settings[f"{genome}_genome_fasta"]
        if config_entry != '__default-ensemble__':
            outfname = config_entry
        elif not cache_path:
            raise CacheUnavailableError('No cache path defined, but default ensembl fasta files should be used.')
        else:
            base_args = (cache_path, 'fasta', 'homo_sapiens', f'{ENSEMBL_RELEASE}_{{genome}}')
            outfname = os.path.join(*base_args, 'Homo_sapiens.{genome}.dna.primary_assembly.fa.gz')
            if genome == 'hg38':
                wildcards['genome'] = 'GRCh38'
            else:
                wildcards['genome'] = 'GRCh37'
    elif type == 'gtf':
        config_entry = global_settings[f"{genome}_gtf_file"]
        if config_entry != '__default-gencode__':
            outfname = config_entry
        elif not cache_path:
            raise CacheUnavailableError('No cache path defined, but VEP fasta files should be used.')
        else:
            outfname = os.path.join(cache_path, 'static-data', 'gencode.{genome}.v45.gtf.gz')
    elif type == 'genome_info':
        config_entry = global_settings[f'{genome}_genomeInfo_file']
        if config_entry != '__default-UCSC__':
            outfname = config_entry
        elif not cache_path:
            raise CacheUnavailableError('No cache path defined, but VEP fasta files should be used.')
        else:
            outfname = os.path.join(cache_path, 'static-data', 'UCSC_{genome}_chromosome-info.tsv')
    elif type == 'mehari_txdb':
        config_entry = global_settings[f'{genome}_mehari_transcript_db']
        wildcards['genome'] = 'GRCh38' if genome == 'hg38' else 'GRCh37'
        wildcards['MEHARI_DB_VERSION'] = MEHARI_DB_VERSION
        if config_entry != '__cache-default__':
            outfname = config_entry
        elif not cache_path:
            raise CacheUnavailableError('No StemCNV-check cache defined, but mehari-db path is set to use cache path.')
        else:
            outfname = os.path.join(
                cache_path,
                'mehari-db',
                "mehari-data-txs-{genome}-ensembl-and-refseq-{MEHARI_DB_VERSION}.bin.zst"
            )
    elif type == 'dosage_scores':
        config_entry = global_settings['dosage_sensitivity_scores']
        if config_entry != '__cache-default__':
            outfname = config_entry
        elif not cache_path:
            raise CacheUnavailableError('No StemCNV-check cache defined, but dosage_sensivity_data is set to use cache path.')
        else:
            outfname = os.path.join(
                cache_path,
                "Collins_rCNV_2022.dosage_sensitivity_scores.tsv.gz"
            )
    else:
        raise ValueError('Unknown global file type: ' + type)

    if fill_wildcards:
        return outfname.format(**wildcards)
    else:
        return outfname


def get_array_file(filekey, array_name, config, cache_path):
    """Get file path from array definition config, using default cache paths where needed."""

    if filekey not in ('penncnv_pfb_file', 'penncnv_GCmodel_file', 'array_density_file', 'array_gaps_file'):
        raise ValueError('Unknown array file key: ' + filekey)

    genome_build = config['array_definition'][array_name]['genome_version']
    genome_build = 'hg38' if genome_build in ('hg38', 'GRCh38') else 'hg19'

    config_entry = config['array_definition'][array_name][filekey]

    if config_entry == '__cache-default__':
        if not cache_path:
            raise CacheUnavailableError('No cache path defined, but default cache array definition files should be used.')
        base_args = (cache_path, 'array_definitions', array_name)
        if filekey == 'penncnv_pfb_file':
            outfname = os.path.join(*base_args, f'PennCNV-PFB_{genome_build}.pfb')
        elif filekey == 'penncnv_GCmodel_file':
            outfname = os.path.join(*base_args, f'PennCNV-GCmodel_{genome_build}.gcmodel')
        elif filekey == 'array_density_file':
            outfname = os.path.join(*base_args, f'density_{genome_build}.bed')
        elif filekey == 'array_gaps_file':
            outfname = os.path.join(*base_args, f'gaps_{genome_build}.bed')
    else:
        outfname = config_entry

    return outfname
