# -*- coding: utf-8 -*-
"""Helper functions for pipeline"""

def read_sample_table(filename):
  samples = {}
  cols = ('SampleName', 'Chip_Name', 'Sample_ID', 'Sex', 'ReferenceSample')
  with open(filename, 'r') as f:
    header_index = {h: cols.index(h) for h in next(f).split() if h in cols}
    if len(header_index) < 5:
      raise Error('Not all require sample_table columns found: SampleName, Chip_Name, Sample_ID, Sex, ReferenceSample')
    for line in f:
      if line.startswith('#'):
        continue
      line = line.rstrip('\n').split('\t')
      line = {col: line[header_index[col]] for col in cols}
      samples[line['SampleName']] = (line['Chip_Name'],
                                     line['Sample_ID'],
                                     line['Sex'],
                                     line['ReferenceSample']
                                    )
  return samples
      
