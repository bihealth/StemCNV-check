# -*- coding: utf-8 -*-
"""Helper functions for pipeline"""
from py_exceptions import *


def read_sample_table(filename):
    samples = {}
    cols = ('Sample_Name', 'Chip_Name', 'Chip_Pos', 'Sample_ID', 'Sex', 'Reference_Sample')
    with open(filename, 'r') as f:
        header = next(f).rstrip('\n')
        while header.startswith('#'):
            header = next(f).rstrip('\n')
        header = header.split('\t')
        header_index = {col: header.index(col) for col in cols if col in header}
        if len(header_index) < len(cols):
            raise SampletableDefaultColError('Not all required sample_table columns found: ' + ', '.join(cols))
        for line in f:
            if line.startswith('#'):
                continue
            line = line.rstrip('\n').split('\t')
            line = {col: line[header_index[col]] for col in cols}
            samples[line['Sample_Name']] = (line['Chip_Name'],
                                            line['Chip_Pos'],
                                            line['Sample_ID'],
                                            line['Sex'],
                                            line['Reference_Sample']
                                            )
    return samples
