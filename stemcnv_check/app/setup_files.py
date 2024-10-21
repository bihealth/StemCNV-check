import importlib.resources
import os
import shutil
from loguru import logger as logging
from .. import STEM_CNV_CHECK
from ..helpers import read_sample_table

def setup_control_files(args):

    logging.info(f'Creating empty config and sample table files. Config details: {args.config_details}')

    if os.path.exists(args.sample_table) and not args.overwrite:
        logging.info(f"Sample table already exists: {args.sample_table}. Use --overwrite to replace it.")
    elif args.sampletable_format == 'tsv':
        shutil.copyfile(importlib.resources.files(STEM_CNV_CHECK).joinpath('control_files', 'sample_table_example.tsv'), args.sample_table)
    elif args.sampletable_format == 'xlsx':
        import pandas as pd
        source_file = importlib.resources.files(STEM_CNV_CHECK).joinpath('control_files', 'sample_table_example.tsv')
        comment_lines = []
        with open(source_file, 'r') as fin:
            for line in fin:
                if line.startswith('#'):
                    comment_lines.append(line.rstrip())
                else:
                    break
        logging.info('writing sample table: ' + args.sample_table)
        with pd.ExcelWriter(args.sample_table) as writer:
            pd.DataFrame(comment_lines).to_excel(writer, index=False, header=False)
            read_sample_table(source_file, return_type='dataframe').to_excel(writer, index=False)

    if os.path.exists(args.config) and not args.overwrite:
        logging.info(f"Config file already exists: {args.config}. Use --overwrite to replace it.")
    else:
        with importlib.resources.files(STEM_CNV_CHECK).joinpath('control_files', 'default_config.yaml').open() as fin, \
                open(args.config, 'w') as fout:
            modes = ('minimal', 'medium', 'advanced', 'complete')
            use_mode = args.config_details
            write_lines = False
            for line in fin:
                if line.startswith('##!'):
                    write_lines = modes.index(line[3:-1]) <= modes.index(use_mode)
                    continue
                if write_lines:
                    fout.write(line)
