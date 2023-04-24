# -*- coding: utf-8 -*-
"""Helper functions for pipeline"""

#TODO:
# allow commented lines above header?
def read_sample_table(filename):
  samples = {}
  cols = ('SampleName', 'Chip_Name', 'Sample_ID', 'Sex', 'ReferenceSample')
  with open(filename, 'r') as f:
    header = next(f).rstrip('\n').split('\t')
    header_index = {col: header.index(col) for col in cols if col in header}
    if len(header_index) < 5:
      raise Exception('Not all required sample_table columns found: SampleName, Chip_Name, Sample_ID, Sex, ReferenceSample')
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
      
