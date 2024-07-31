import importlib.resources
import os

import ruamel.yaml as ruamel_yaml

from snakemake.cli import main
from loguru import logger as logging
from .. import STEM_CNV_CHECK
from ..helpers import make_PennCNV_sexfile, make_apptainer_args

def run_stemcnv_check_workflow(args):
    # Ensure that sexfile for PennCNV exists
    make_PennCNV_sexfile(args)
    with open(args.config) as f:
        yaml = ruamel_yaml.YAML(typ='safe')
        config = yaml.load(f)

    argv = [
        "-s", str(importlib.resources.files(STEM_CNV_CHECK).joinpath('rules', 'StemCNV-check.smk')),
        "--printshellcmds", "--rerun-incomplete",
        "--sdm", "conda", "apptainer",
        "--apptainer-args", make_apptainer_args(config),
        #"--conda-frontend", args.conda_frontend,
        ]
    if args.directory:
        argv += ["--directory", args.directory]

    basedir = args.directory if args.directory else os.getcwd()
    argv += [
        '--configfile', str(importlib.resources.files(STEM_CNV_CHECK).joinpath('control_files', 'default_config.yaml')),
        args.config,
        '--config', f'sample_table={args.sample_table}',
        f'basedir={basedir}',
        f'configfile={args.config}',
        f'target={args.target}',
    ]
    #FIXME: use a clearer local vs cluster submission
    if args.cluster_profile:
        argv += [
            "--profile", args.cluster_profile,
            "-j", str(args.jobs)
        ]
    else:
        argv += [
            '--cores', str(args.local_cores)
        ]

    if args.snake_options:
        argv += args.snake_options

    return main(argv)
