import importlib.resources
import ruamel.yaml as ruamel_yaml

from snakemake import main
from loguru import logger as logging
from .. import STEM_CNV_CHECK
from ..helpers import make_PennCNV_sexfile, make_singularity_args

def run_stemcnv_check_workflow(args):
    # Ensure that sexfile for PennCNV exists
    make_PennCNV_sexfile(args)
    with open(args.config) as f:
        yaml = ruamel_yaml.YAML(typ='safe')
        config = yaml.load(f)

    argv = [
        "-s", importlib.resources.files(STEM_CNV_CHECK).joinpath('rules', 'StemCNV-check.smk'),
        "-p", "--rerun-incomplete",
        "--use-conda", "--conda-frontend", args.conda_frontend,
    ]
    if args.no_singularity:
        logging.warning(
            "Running without singularity containers is a legacy feature, pipeline will fail without local PennCNV installation based on the install.sh script!")
        use_singularity = False
    else:
        argv += ["--use-singularity",
                 "--singularity-args", make_singularity_args(config)
                 ]
        use_singularity = True

    argv += [
        '-d', args.directory,
        '--configfile', importlib.resources.files(STEM_CNV_CHECK).joinpath('control_files', 'default_config.yaml'),
        '--config', f'sample_table={args.sample_table}',
        #f'snakedir={SNAKEDIR}',
        f'basedir={args.directory}',
        f'configfile={args.config}',
        f'target={args.target}',
        f'use_singularity={use_singularity}'
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
