import importlib.resources
import os
import re
import shutil
import ruamel.yaml as ruamel_yaml
from loguru import logger as logging
from stemcnv_check import STEM_CNV_CHECK
from stemcnv_check.helpers import read_sample_table

def setup_control_files(args):

    logging.info(f'Creating empty config and sample table files. Config details: {args.config_details}')

    if os.path.exists(args.sample_table) and not args.overwrite:
        logging.warning(f"Sample table already exists: {args.sample_table}. Use --overwrite to replace it.")
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
        with pd.ExcelWriter(args.sample_table) as writer:
            pd.DataFrame(comment_lines).to_excel(writer, index=False, header=False)
            read_sample_table(source_file).to_excel(writer, index=False)

    if os.path.exists(args.config) and not args.overwrite:
        logging.warning(f"Config file already exists: {args.config}. Use --overwrite to replace it.")
    else:
        # Load defined labels to autofill comments in config file
        yaml = ruamel_yaml.YAML(typ='safe')
        with importlib.resources.files(STEM_CNV_CHECK).joinpath('control_files').joinpath(
                'label_name_definitions.yaml').open() as f:
            defined_label_names = yaml.load(f)
        comment_format = {
            k: ', '.join(v) for k, v in defined_label_names.items() if isinstance(v, list)
        }
        comment_format.update({
            k+'_names': ', '.join(v.keys()) for k, v in defined_label_names.items() if isinstance(v, dict)
        })
        comment_format.update({
            k+'_values': ', '.join(v.values()) for k, v in defined_label_names.items() if isinstance(v, dict)
        })

        with importlib.resources.files(STEM_CNV_CHECK).joinpath('control_files', 'default_config.yaml').open() as fin, \
                open(args.config, 'w') as fout:
            modes = ('minimal', 'medium', 'advanced', 'complete')
            use_mode = args.config_details
            write_lines = False
            for line in fin:
                # 'Uncomment' lines starting with '#!' in default
                if line.startswith('#!'):
                    line = line[2:]
                # modify lines with '###!' to strip characters from the end
                if '###!strip:' in line:
                    n_strip = int(line[line.index('###!strip:')+10:])
                    line = line[:line.index('###!strip:')-n_strip] + '\n'
                # Format normal comment lines as strings
                if re.match(r'\s+#', line):
                    line = line.format(**comment_format)
                # '##!' lines are used for control of written out detail
                if line.startswith('##!'):
                    write_lines = modes.index(line[3:-1]) <= modes.index(use_mode)
                    continue
                if write_lines:
                    fout.write(line)
