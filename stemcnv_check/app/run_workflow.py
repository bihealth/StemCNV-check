import importlib.resources
import os
import ruamel.yaml as ruamel_yaml

from snakemake.cli import main
from loguru import logger as logging
from .. import STEM_CNV_CHECK
from ..helpers import load_config, make_apptainer_args, get_cache_dir, get_vep_cache_path

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
            "--conda-prefix", cache_path
        ]
    # Some disc space could be cleared by using these options:
    # --cleanup-containers
    # --conda-cleanup-envs

    # Define / overwrite place-holder values for VEP downloaded data
    vep_cache_path = get_vep_cache_path(config['settings']['VEP_annotation']['VEP_cache_path'], cache_path)
    hg19_fasta = config['global_settings']['hg19_genome_fasta'].replace(
        '__use-vep__',
        os.path.join(vep_cache_path, 'fasta', 'homo_sapiens', '112_GRCh37', 'Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz')
    )
    hg38_fasta = config['global_settings']['hg38_genome_fasta'].replace(
        '__use-vep__',
        os.path.join(vep_cache_path, 'fasta', 'homo_sapiens', '112_GRCh38', 'Homo_sapiens.GRCh38.dna.toplevel.fa.gz')
    )

    basedir = args.directory if args.directory else os.getcwd()
    argv += [
        '--configfile', str(importlib.resources.files(STEM_CNV_CHECK).joinpath('control_files', 'default_config.yaml')),
        args.config,
        '--config', f'sample_table={args.sample_table}',
        f'basedir={basedir}',
        f'configfile={args.config}',
        f'target={args.target}',
        f'use_vep_cache={vep_cache_path}',
        f"global_settings='{{hg19_genome_fasta: {hg19_fasta}, hg38_genome_fasta: {hg38_fasta}}}'"
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
