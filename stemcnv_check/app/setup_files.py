import importlib.resources
import os
import shutil
from loguru import logger as logging

from .. import STEM_CNV_CHECK

def setup_control_files(args):

    logging.info(f'Creating empty config and sample table files. Config details: {args.config_details}')

    if os.path.exists(args.sample_table) and not args.overwrite:
        logging.info(f"Sample table already exists: {args.sample_table}. Use --overwrite to replace it.")
    else:
        shutil.copyfile(importlib.resources.files(STEM_CNV_CHECK).joinpath('control_files', 'sample_table_example.tsv'), args.sample_table)

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
