import importlib.resources
import ruamel.yaml as ruamel_yaml
import os
from snakemake import main
from loguru import logger as logging
from .. import STEM_CNV_CHECK
from ..helpers import make_PennCNV_sexfile, make_singularity_args, get_cache_dir, get_vep_cache_path, config_extract

def run_stemcnv_check_workflow(args):
    # Ensure that sexfile for PennCNV exists
    make_PennCNV_sexfile(args)
    with open(args.config) as f:
        yaml = ruamel_yaml.YAML(typ='safe')
        config = yaml.load(f)
    default_config = importlib.resources.files(STEM_CNV_CHECK).joinpath('control_files', 'default_config.yaml')
    with default_config.open() as f:
        DEF_CONFIG = yaml.load(f)

    argv = [
        "-s", str(importlib.resources.files(STEM_CNV_CHECK).joinpath('rules', 'StemCNV-check.smk')),
        "-p", "--rerun-incomplete",
        "--use-conda", "--conda-frontend", args.conda_frontend,
    ]

    # make snakemake write conda & singularity files to an "external" cache-path, NOT individual project paths
    # This saves disk space when using multiple projects with the same conda envs
    cache_path = get_cache_dir(args)
    if cache_path:
        # snakemake.main can not deal with Path objects
        cache_path = str(cache_path)
        argv += [
            "--conda-prefix", cache_path
        ]
    # Some disc space could be cleared by using these options:
    # --cleanup-containers
    # --conda-cleanup-envs
    # Use this to overwrite the place-holder in the given config
    vep_cache_path = get_vep_cache_path(
        config_extract(('static_data', 'VEP_cache_path'), config, DEF_CONFIG),
        cache_path
    )

    if args.no_singularity:
        logging.warning(
            "Running without singularity containers is a legacy feature, pipeline will fail without local PennCNV installation based on the install.sh script!")
        use_singularity = False
    else:
        argv += ["--use-singularity",
                 "--singularity-args", make_singularity_args(config)
                 ]
        use_singularity = True
        if cache_path:
            argv += [
                "--singularity-prefix", cache_path
            ]

    if args.directory:
        argv += ["--directory", args.directory]


    argv += [
        '--configfile', str(importlib.resources.files(STEM_CNV_CHECK).joinpath('control_files', 'default_config.yaml')),
        args.config,
        '--config', f'sample_table={args.sample_table}',
        # f'basedir={args.directory}',
        f'configfile={args.config}',
        f'target={args.target}',
        f'use_singularity={use_singularity}',
        f'genome_build={args.genome}',
        f'static_data={{"VEP_cache_path":{vep_cache_path}}}'
    ]

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
