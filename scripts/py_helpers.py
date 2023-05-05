# -*- coding: utf-8 -*-
"""Helper functions for pipeline"""
from .py_exceptions import *


def read_sample_table(filename):
    samples = []
    cols = ('Sample_ID', 'Chip_Name', 'Chip_Pos', 'Sex', 'Reference_Sample')
    with open(filename, 'r') as f:
        header = next(f).rstrip('\n')
        while header.startswith('#'):
            header = next(f).rstrip('\n')
        header = header.split('\t')
        header_index = {col: header.index(col) for col in cols if col in header}
        if len(header_index) < len(cols):
            missing = [c for c in cols if c not in header_index.keys()]
            raise SampletableDefaultColError('Not all required sample_table columns found. Missing columns: ' + ', '.join(missing))
        for line in f:
            line = line.rstrip('\n')
            if not line:
                continue
            if line.startswith('#'):
                continue
            line = line.split('\t')
            line = {col: line[header_index[col]] for col in cols}
            samples.append([line['Sample_ID'],
                            line['Chip_Name'],
                            line['Chip_Pos'],
                            line['Sex'],
                            line['Reference_Sample']
                            ])
    return samples
