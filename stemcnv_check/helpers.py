# -*- coding: utf-8 -*-
"""Helper functions for pipeline"""
import importlib.resources
import os
import ruamel.yaml as ruamel_yaml
from pathlib import Path
from pydantic.v1.utils import deep_update
from . import STEM_CNV_CHECK, mehari_db_version
from .exceptions import SampleConstraintError, ConfigValueError
from collections import OrderedDict
from loguru import logger as logging


def read_sample_table(filename, with_opt=False):
    samples = []
    cols = ['Sample_ID', 'Chip_Name', 'Chip_Pos', 'Sex', 'Reference_Sample']
    # req_cols = cols[:]
    # if with_opt:
    #     cols += ['Sample_Group', 'Regions_of_Interest']
    with open(filename, 'r') as f:
        header = next(f).rstrip('\n')
        while header.startswith('#'):
            header = next(f).rstrip('\n')
        header = header.split('\t')
        if not all(col in header for col in cols):
            missing = [c for c in cols if c not in header]
            raise SampleConstraintError('Not all required sample_table columns found. Missing columns: ' + ', '.join(missing))
        # add optinal cols
        if with_opt:
            cols = cols + [c for c in header if c not in cols]
        for line in f:
            line = line.rstrip('\n')
            if not line:
                continue
            if line.startswith('#'):
                continue
            line = line.split('\t')
            # reorder line
            try:
                line_ordered = [line[header.index(val)] for val in cols]
            except IndexError:
                logging.error("Could not extract proper columns from this line:\n" + '\t'.join(line))
                raise
            if with_opt:
                samples.append(OrderedDict((name, val) for name, val in zip(cols, line_ordered)))
            else:
                samples.append(line_ordered)

    return samples


def make_apptainer_args(config, tmpdir=None, not_existing_ok=False):
    """Collect all outside filepaths that need to be bound inside container"""

    bind_points = [
        (config['data_path'], '/outside/data'),
        (config['raw_data_folder'], '/outside/rawdata'),
        (config['log_path'], '/outside/logs'),
        (importlib.resources.files(STEM_CNV_CHECK), '/outside/snakedir'),
    ]
    if tmpdir is not None:
        bind_points.append((tmpdir, '/outside/tmp'))

    for name, file in config['static_data'].items():
        # Can only mount existing files
        if not os.path.isfile(file):
            if not_existing_ok or not file:
                continue
            else:
                raise FileNotFoundError(f"Static data file '{file}' does not exist.")
        bind_points.append((file, '/outside/static/{}'.format(os.path.basename(file))))

    bind_point_str = "-B " + ','.join(f"{host}:{cont}" for host, cont in bind_points)
    logging.debug("Binding points for apptainer: " + str(bind_point_str))

    return bind_point_str


def collect_SNP_cluster_ids(sample_id, config_extra_samples, sample_data_full):
    ids = []
    # '__[column]' entries: take all sample_ids with the same value in '[column]'
    col_val_match = [sampledef[2:] for sampledef in config_extra_samples if sampledef[:2] == '__']
    for col in col_val_match:
        if col not in sample_data_full[0].keys():
            raise ConfigValueError('Config for SNP clustering refers to non-existing column: ' + col)
        match_val = [dictline[col] for dictline in sample_data_full if dictline['Sample_ID'] == sample_id][0]
        ids += [dictline['Sample_ID'] for dictline in sample_data_full if dictline[col] == match_val]
    # '_[column]' entry: take all sample_ids from '[column]'
    id_cols = [sampledef[1:] for sampledef in config_extra_samples if sampledef[0] == '_' and sampledef[:2] != '__']
    for col in id_cols:
        if col not in sample_data_full[0].keys():
            raise ConfigValueError('Config for SNP clustering refers to non-existing column: ' + col)
        ids += [dictline[col] for dictline in sample_data_full if dictline['Sample_ID'] == sample_id][0].split(',')
    # other entries: assume they are sample_ids & use them as is
    ids += [sampledef for sampledef in config_extra_samples if sampledef[0] != '_']

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
        if path.is_dir() and os.access(path, os.W_OK):
            logging.debug(f"Using existing cache directory: {path}")
            return path
        else:
            try:
                path.mkdir()
                logging.debug(f"Created cache directory: {path}")
                return path
            except PermissionError:
                logging.debug(f"Could not create cache in {path}")
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
            logging.warning(f"Could not use cache directory '{args.cache_path}'")
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

    # Nothing worked, return None & don't use specific cache
    logging.info("Cache directory is not usable! Conda and docker images will be stored in snakemake project directory")
    return None


def get_mehari_db_file(config_entry, cache_path, genome_version):
    # If given use specific path
    if config_entry != '__cache-default__':
        return config_entry
    # If cache should be used but is missing, throw error
    elif not cache_path:
        raise ConfigValueError('No StemCNV-check cache defined, but mehari-db path is set to use cache path.')
    # Use same location as cache for conda & docker from workflow
    else:
        genome_version = 'GRCh38' if genome_version in ("hg38", "GRCh38") else 'GRCh37'
        return os.path.join(
            cache_path,
            'mehari-db',
            f"mehari-data-txs-{genome_version}-ensembl-{mehari_db_version}.bin.zst"
        )

