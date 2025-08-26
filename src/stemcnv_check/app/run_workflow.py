import importlib.resources
import os
import ruamel.yaml as ruamel_yaml

from snakemake.cli import main
from loguru import logger as logging
from stemcnv_check import STEM_CNV_CHECK
from stemcnv_check.helpers import load_config, make_apptainer_args, get_cache_dir, get_cache_array_definition

def run_stemcnv_check_workflow(args, is_wsl):

    config = load_config(args)
    cache_path = get_cache_dir(args, config)

    argv = [
        "-s", str(importlib.resources.files(STEM_CNV_CHECK).joinpath('rules', 'StemCNV-check.smk')),
        "--printshellcmds", "--rerun-incomplete",
        "--sdm", "conda", "apptainer",
        "--apptainer-args", make_apptainer_args(config, cache_path, extra_bind_args=args.bind_points),
        ]
    if args.directory:
        argv += ["--directory", args.directory]

    # make snakemake write conda & singularity files to an "external" cache-path, NOT individual project paths
    # This saves disk space when using multiple projects with the same conda envs
    if cache_path:
        # snakemake.main can not deal with Path objects
        cache_path = str(cache_path)
        argv += [
            "--conda-prefix", cache_path,
            "--apptainer-prefix", cache_path
        ]

    basedir = args.directory if args.directory else os.getcwd()
    argv += [
        '--configfile', str(importlib.resources.files(STEM_CNV_CHECK).joinpath('control_files', 'default_config.yaml')),
    ]
    # Only add global config if it exists (inbetween defaults and user config)
    # Adding an emtpy string will cause issues for snakemake
    global_config = get_cache_array_definition(cache_path)
    if global_config and os.path.isfile(global_config):
        argv += [global_config]
    argv += [
        args.config,
        '--config',
        f'sample_table={args.sample_table}',
        f"column_remove_regex={args.column_remove_regex}",
        f'basedir={basedir}',
        f'configfile={args.config}',
        f'target={args.target}',
        f'cache_path={cache_path}',
        f'verbose_level={args.verbose}',
        f'is_wsl={is_wsl}',
    ]
    if args.collate_date:
        argv += [f'collate_date={args.collate_date}']

    argv += [
        '--cores', str(args.local_cores)
    ]
    if args.memory_mb:
        argv += ["--resources", f"mem_mb={args.memory_mb}"]

    if args.snake_options:
        argv += args.snake_options

    logging.debug(f"Running snakemake with this command:\nsnakemake{' '.join(argv)}")
    return main(argv)
