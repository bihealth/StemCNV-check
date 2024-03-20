# -*- coding: utf-8 -*-
"""Helper functions for pipeline"""
from .py_exceptions import *
from collections import OrderedDict
from collections.abc import MutableMapping
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
            raise SampletableDefaultColError('Not all required sample_table columns found. Missing columns: ' + ', '.join(missing))
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
            line_ordered = [line[header.index(val)] for val in cols]
            if with_opt:
                samples.append(OrderedDict((name, val) for name, val in zip(cols, line_ordered)))
            else:
                samples.append(line_ordered)

    return samples


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


def config_extract(entry_kws, config, def_config, verbose=False):
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
            if verbose:
                print('Warning: Using config default values for: ' + ' : '.join(used_entries))
        else:
            logging.warning(' : '.join(used_entries) + " is not a valid config entry or has been deprecated", ConfigKeyWarning)
            return None

    return subconfig


def flatten(dictionary, parent_key='', separator=':'):
    items = []
    for key, value in dictionary.items():
        new_key = parent_key + separator + key if parent_key else key
        if isinstance(value, MutableMapping):
            items.extend(flatten(value, new_key, separator=separator).items())
        else:
            items.append((new_key, value))
    return dict(items)
