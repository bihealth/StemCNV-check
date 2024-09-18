import importlib.resources
import os
import ruamel.yaml as ruamel_yaml

from snakemake.cli import main
from loguru import logger as logging
from .. import STEM_CNV_CHECK
from ..helpers import load_config, make_apptainer_args, get_cache_dir

def run_stemcnv_check_workflow(args):

    config = load_config(args.config)

    argv = [
        "-s", str(importlib.resources.files(STEM_CNV_CHECK).joinpath('rules', 'StemCNV-check.smk')),
        "--printshellcmds", "--rerun-incomplete",
        "--sdm", "conda", "apptainer",
        "--apptainer-args", make_apptainer_args(config),
        #"--conda-frontend", args.conda_frontend,
        ]
    if args.directory:
        argv += ["--directory", args.directory]

    # make snakemake write conda & singularity files to an "external" cache-path, NOT individual project paths
    # This saves disk space when using multiple projects with the same conda envs
    cache_path = get_cache_dir(args, config)
    if cache_path:
        # snakemake.main can not deal with Path objects
        cache_path = str(cache_path)
        argv += [
            "--conda-prefix", cache_path,
            "--apptainer-prefix", cache_path
        ]
    # Some disc space could be cleared by using these options:
    # --cleanup-containers
    # --conda-cleanup-envs


    basedir = args.directory if args.directory else os.getcwd()
    argv += [
        '--configfile', str(importlib.resources.files(STEM_CNV_CHECK).joinpath('control_files', 'default_config.yaml')),
        args.config,
        '--config',
        f'sample_table={args.sample_table}',
        f'basedir={basedir}',
        f'configfile={args.config}',
        f'target={args.target}',
        f'cache_path={cache_path}',
        f'verbose_level={args.verbose}',
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
