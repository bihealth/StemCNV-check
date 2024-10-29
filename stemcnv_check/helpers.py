# -*- coding: utf-8 -*-
"""Helper functions for pipeline"""
import importlib.resources
import pandas as pd
import os
import re
import ruamel.yaml as ruamel_yaml
from pathlib import Path
from pydantic.v1.utils import deep_update
from . import STEM_CNV_CHECK, mehari_db_version, VEP_version
from .exceptions import SampleConstraintError, ConfigValueError, CacheUnavailableError
from collections import OrderedDict
from loguru import logger as logging


def read_sample_table(filename, name_remove_regex=None, return_type='list'):
    """Read sample table from file, return either a list of the required col entires,
    a list with Ordereddict of all rows (colnames as keys) or a pandas dataframe"""
    
    if return_type not in ('list', 'list_withopt', 'dataframe'):
        raise ValueError('Invalid return_type: ' + return_type)

    if str(filename).endswith('.xlsx'):
        def reader_func(f): return pd.read_excel(f, comment='#')
    elif str(filename).endswith('.tsv') or str(filename).endswith('.txt'):
        def reader_func(f): return pd.read_csv(f, sep='\t', comment='#')
    elif str(filename).endswith('.csv'):
        def reader_func(f): return pd.read_csv(f, comment='#')
    else:
        raise ValueError('Unknown file format for sample table: ' + filename)

    req_cols = ['Sample_ID', 'Chip_Name', 'Chip_Pos', 'Array_Name', 'Sex', 'Reference_Sample']

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


def make_apptainer_args(config, cache_path, tmpdir=None, not_existing_ok=False):
    """Collect all outside filepaths that need to be bound inside container"""

    bind_points = [
        (config['data_path'], '/outside/data'),
        (config['raw_data_folder'], '/outside/rawdata'),
        (config['log_path'], '/outside/logs'),
        (importlib.resources.files(STEM_CNV_CHECK), '/outside/snakedir'),
    ]

    # When apptainer is used these should always be present
    used_genomes = set(array['genome_version'] for name, array in config['array_definition'].items() if name != '_default_')
    for global_file in ('fasta', 'gtf', 'genome_info', 'mehari_txdb'):
        for genome_version in used_genomes:
            static_file = get_global_file(global_file, genome_version, config['global_settings'], cache_path)
            # if not os.path.isfile(static_file):
            #     if not_existing_ok or not static_file:
            #         continue
            #     else:
            #         raise FileNotFoundError(f"Static data file '{static_file}' does not exist.")
            bind_points.append((static_file, '/outside/static/' + os.path.basename(static_file)))

    if tmpdir is not None:
        bind_points.append((tmpdir, '/outside/tmp'))

    for array in config['array_definition'].keys():
        if array == '_default_':
            continue
        for name, file in config['array_definition'][array].items():
            if name == 'genome_version':
                continue
            # Can only mount existing files
            if not os.path.isfile(file):
                if not_existing_ok or not file:
                    continue
                else:
                    raise FileNotFoundError(f"Static data file '{file}' does not exist.")
            bind_points.append((file, '/outside/{array}/{fname}'.format(array=array, fname=os.path.basename(file))))

    # Sort bind_points to make testing easier (loops over config_dict add them in a non-deterministic order)
    bind_point_str = "-B " + ','.join(f"'{host}':'{cont}'" for host, cont in sorted(bind_points, key=lambda x: x[1]))
    logging.debug("Binding points for apptainer: " + str(bind_point_str))

    return bind_point_str


def collect_SNP_cluster_ids(sample_id, config_extra_samples, sample_data_df):
    ids = set()
    # '__[column]' entries: take all sample_ids with the same value in '[column]'
    col_val_match = [sampledef[2:] for sampledef in config_extra_samples if sampledef[:2] == '__']
    for col in col_val_match:
        if col not in sample_data_df.columns:
            raise ConfigValueError('Config for SNP clustering refers to non-existing column: ' + col)
        match_val = sample_data_df[col].loc[sample_id]
        ids.update(sample_data_df.set_index(col)['Sample_ID'].loc[match_val])
    # '_[column]' entry: take all sample_ids from '[column]'
    id_cols = [sampledef[1:] for sampledef in config_extra_samples if sampledef[0] == '_' and sampledef[:2] != '__']
    for col in id_cols:
        if col not in sample_data_df.columns:
            raise ConfigValueError('Config for SNP clustering refers to non-existing column: ' + col)
        ids.update(sample_data_df[col].loc[sample_id].split(','))
    # other entries: assume they are sample_ids & use them as is
    ids.update([sampledef for sampledef in config_extra_samples if sampledef[0] != '_'])
    # remove the original sample_id
    if sample_id in ids:
        ids.remove(sample_id)

    # Check that all samples belong to the same array,
    # remove all that don't match
    sample_array = sample_data_df['Array_Name'].loc[sample_id]
    wrong_array = set()
    for id in ids:
        if sample_data_df['Array_Name'].loc[id] != sample_array:
            wrong_array.add(id)
    if wrong_array:
        logging.warning(f"Samples {', '.join(wrong_array)} do not belong to the same array as {sample_id}. They will be excluded from the SNP clustering.")
        ids -= wrong_array

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


def load_config(configfile, defaults=True):
    yaml = ruamel_yaml.YAML(typ='safe')
    with open(configfile) as f:
        config = yaml.load(f)

    if defaults:
        config_def_file = importlib.resources.files(STEM_CNV_CHECK).joinpath('control_files', 'default_config.yaml')
        with config_def_file.open() as f:
            default_config = yaml.load(f)
        return deep_update(default_config, config)
    else:
        return config


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
        auto_paths = [
            ('config-defined', Path(config['global_settings']['cache_dir']).expanduser()),
            #('install-dir', importlib.resources.files(STEM_CNV_CHECK).joinpath('.stemcnv-check-cache')),
            ('home', Path('~/.cache/stemcnv-check').expanduser())
        ]
        for name, path in auto_paths:
            cache_path = check_path_usable(path)
            if cache_path:
                return cache_path
            else:
                logging.debug(f"Could not use cache directory '{path}'")

    # Nothing worked, return None & don't use specific cache
    logging.warning("Cache directory is not usable! Trying to run without cache. Conda and docker images will be stored in snakemake project directory")
    return None


def get_global_file(type, genome_version, global_settings, cache_path, fill_wildcards = True):
    """Get path to files defined in the global settings. All of these have an internally defined default.
    Supported types are: fasta, gtf, genome_info, mehari_txdb"""

    genome = 'hg38' if genome_version in ("hg38", "GRCh38") else 'hg19'
    wildcards = {'genome': genome_version}

    if type == 'fasta':
        config_entry = global_settings[f"{genome}_genome_fasta"]
        if config_entry != '__use-vep__':
            outfname = config_entry
        elif not cache_path:
            raise CacheUnavailableError('No cache path defined, but VEP fasta files should be used.')
        else:
            base_args = (cache_path, 'fasta', 'homo_sapiens', f'{VEP_version}_{{genome}}')
            if genome == 'hg38':
                wildcards['genome'] = 'GRCh38'
                outfname = os.path.join(*base_args, 'Homo_sapiens.{genome}.dna.toplevel.fa.gz')
            else:
                wildcards['genome'] = 'GRCh37'
                outfname = os.path.join(*base_args, 'Homo_sapiens.{genome}.75.dna.primary_assembly.fa.gz')

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
        config_entry = global_settings['mehari_transcript_db']
        wildcards['genome'] = 'GRCh38' if genome == 'hg38' else 'GRCh37'
        wildcards['mehari_db_version'] = mehari_db_version
        if config_entry != '__cache-default__':
            outfname = config_entry
        elif not cache_path:
            raise CacheUnavailableError('No StemCNV-check cache defined, but mehari-db path is set to use cache path.')
        else:
            outfname = os.path.join(
                cache_path,
                'mehari-db',
                "mehari-data-txs-{genome}-ensembl-{mehari_db_version}.bin.zst"
            )
    else:
        raise ValueError('Unknown global file type: ' + type)

    if fill_wildcards:
        return outfname.format(**wildcards)
    else:
        return outfname