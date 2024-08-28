import os
import sys
import tempfile
import ruamel.yaml as ruamel_yaml
from pathlib import Path

from snakemake.api import SnakemakeApi
from snakemake.settings.types import ResourceSettings, ConfigSettings, DeploymentSettings, DAGSettings, OutputSettings, DeploymentMethod
from .check_input import check_config
from .. import STEM_CNV_CHECK, VEP_version
from ..helpers import make_apptainer_args, read_sample_table, get_cache_dir, get_mehari_db_file, load_config
from loguru import logger as logging

import importlib.resources

def create_missing_staticdata(args):
    """
    Use the staticdata_creation.smk workflow to generate missing static files.
    This includes download of fasta and gtf files unless they are already defined in the config.
    Runs first steps of the main StemCNV-check workflow to generate a vcf file if none are present.
    This is done because the vcf files generated with gtc2vcf contain the full information about GT
    (background) counts from the egt clusterfile, that can be used to make the pfb file for PennCNV.
    """
    logging.info('Starting to check for missing static data files ...')
    # bpm & egt files are needed
    check_config(args.config, args.sample_table, required_only=True)
    config = load_config(args.config)
    cache_path = get_cache_dir(args, config)

    array_name = config['array_name']
    genome_build = config['genome_version']
    genome_build = 'hg38' if genome_build in ('hg38', 'GRCh38') else 'hg19'
    ensemble_genome = 'GRCh38' if genome_build == 'hg38' else 'GRCh37'
    mehari_db_file = get_mehari_db_file(config['global_settings']['mehari_transcript_db'], cache_path, ensemble_genome)

    static_snake_config = {
        'TMPDIR': '',
        'genome': genome_build,
        'vcf_input_file': '',
        'density_windows': config['settings']['array_attribute_summary']['density.windows'],
        'min_gap_size': config['settings']['array_attribute_summary']['min.gap.size'],
        'mehari_db_path': os.path.dirname(mehari_db_file),
        'fasta_path': os.path.join(cache_path, 'fasta'),
    }
    # Required static-data files
    static_files = {
        'genome_gtf_file', 'penncnv_pfb_file', 'penncnv_GCmodel_file',
        'genomeInfo_file', 'array_density_file', 'array_gaps_file'
    }
    # Update static_snake_config with existing and target file paths
    existing_files = {file for file in static_files
                      if file in config['static_data'] and os.path.isfile(config['static_data'][file])}
    missing_files = static_files - existing_files
    static_snake_config.update({
        file: config['static_data'][file] for file in existing_files
    })
    static_snake_config.update({
        file: getattr(args, file.lower()).format(genome=genome_build, array=array_name) for file in missing_files
    })

    # Check if mehari db is present
    # FIXME (future): gtf should move here, when array specific files are done as a set
    global_files = set()
    if not os.path.isfile(mehari_db_file):
        os.makedirs(os.path.dirname(mehari_db_file), exist_ok=True)
        global_files.add(mehari_db_file)
    # Check if fasta file is present
    genome_fasta = config['global_settings'][f'{genome_build}_genome_fasta']
    if genome_fasta == '__use-vep__':
        genome_fasta = os.path.join(cache_path, 'fasta',
                                    'homo_sapiens', f'{VEP_version}_{ensemble_genome}',
                                    'Homo_sapiens.GRCh38.dna.toplevel.fa.gz' if ensemble_genome == 'GRCh38' else 
                                    'Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz'
                                    )
    if not os.path.isfile(genome_fasta):
        global_files.add(genome_fasta)

    if not missing_files and not global_files:
        logging.info('All static files are present')
        return 0

    if global_files:
        logging.info(f'Setting up genome fasta file and mehari-db for {genome_build}')
        with (SnakemakeApi(output_settings=OutputSettings(printshellcmds=True)) as api,
              tempfile.TemporaryDirectory() as tmpdir):
            result = (
                api.workflow(
                    snakefile=importlib.resources.files(STEM_CNV_CHECK).joinpath('rules', 'staticdata_creation.smk'),
                    workdir=args.directory,
                    resource_settings=ResourceSettings(
                        cores=args.local_cores,
                        local_cores=args.local_cores,
                        nodes=args.jobs,
                    ),
                    config_settings=ConfigSettings(
                        config=dict(static_snake_config, **{'TMPDIR': tmpdir})
                    ),
                    deployment_settings=DeploymentSettings(
                        deployment_method=DeploymentMethod.parse_choices_set({'conda', 'apptainer'}),
                        apptainer_args=make_apptainer_args(config, not_existing_ok=True),
                        conda_prefix=cache_path,
                        apptainer_prefix=cache_path,
                    ),
                )
                .dag(
                    DAGSettings(targets=global_files,
                                force_incomplete=True)
                )
                .execute_workflow()
            )

    # Check if vcf file is present, generate one if none are
    sample_data = read_sample_table(args.sample_table)
    datapath = config['data_path']
    filter_settings = config['settings']['default-filter-set']
    #FIXME (future): check if annotation is enabled/disabled
    vcf_files = [os.path.join(datapath, f"{sample_id}", f"{sample_id}.annotated-SNP-data.{filter_settings}-filter.vcf.gz") for
                 sample_id, _, _, _, _ in sample_data]
    vcf_present = [vcf for vcf in vcf_files if os.path.exists(vcf)]

    if vcf_present:
        use_vcf = vcf_present[0]
    else:
        use_vcf = vcf_files[0]
        logging.info(f'Running first steps of StemCNV-check to get a vcf file: {use_vcf}')
        with SnakemakeApi(output_settings=OutputSettings(printshellcmds=True)) as api:
            result = (
                    api.workflow(
                        snakefile=importlib.resources.files(STEM_CNV_CHECK).joinpath('rules', 'StemCNV-check.smk'),
                        workdir=args.directory,
                        resource_settings=ResourceSettings(
                            cores=args.local_cores,
                            local_cores=args.local_cores,
                            nodes=args.jobs,
                        ),
                        config_settings=ConfigSettings(
                            configfiles=[
                                importlib.resources.files(STEM_CNV_CHECK).joinpath('control_files', 'default_config.yaml'),
                                Path(args.config),
                            ],
                            config={
                                'sample_table': args.sample_table,
                                'basedir': args.directory,
                                # 'configfile': args.config,
                                # 'target': 'SNP-probe-data',
                                'cache_path': str(cache_path),
                            }
                        ),
                        deployment_settings=DeploymentSettings(
                            deployment_method=DeploymentMethod.parse_choices_set(['conda', 'apptainer']),
                            conda_prefix=cache_path,
                            apptainer_prefix=cache_path,
                            apptainer_args=make_apptainer_args(config, not_existing_ok=True),
                        ),
                    )
                    .dag(
                        DAGSettings(targets={use_vcf},
                                    force_incomplete=True)
                    )
                    .execute_workflow()
                )

    # Run snakemake to generate missing static files
    logging.info(f'Running staticdata creation workflow')

    with (SnakemakeApi(output_settings=OutputSettings(printshellcmds=True)) as api,
          tempfile.TemporaryDirectory() as tmpdir):
        result = (
            api.workflow(
                snakefile=importlib.resources.files(STEM_CNV_CHECK).joinpath('rules', 'staticdata_creation.smk'),
                workdir=args.directory,
                resource_settings=ResourceSettings(
                    cores=args.local_cores,
                    local_cores=args.local_cores,
                    nodes=args.jobs,
                ),
                config_settings=ConfigSettings(
                    config=dict(static_snake_config, **{'TMPDIR': tmpdir, 'vcf_input_file': use_vcf})
                ),
                deployment_settings=DeploymentSettings(
                    deployment_method=DeploymentMethod.parse_choices_set(['conda', 'apptainer']),
                    conda_prefix=cache_path,
                    apptainer_prefix=cache_path,
                    apptainer_args=make_apptainer_args(config, tmpdir=tmpdir, not_existing_ok=True),
                )
            )
            .dag(
                DAGSettings(targets={static_snake_config[file] for file in missing_files},
                            force_incomplete=True)
            )
            .execute_workflow()
        )

    update_str = '\n'.join([f"  {entry}: {file}" for entry, file in static_snake_config.items()
                            if entry in missing_files and
                            entry in config['static_data'] and file != config['static_data'][entry]])

    if args.edit_config_inplace and update_str:
        logging.info('Missing static files generated. Updating config file with new static data entries:\n' + update_str)
        user_config = load_config(args.config, False)
        for file in missing_files:
            user_config['static_data'][file] = static_snake_config[file]
        yaml = ruamel_yaml.YAML()
        with open(args.config, 'w') as f:
            yaml.dump(user_config, f)
    elif update_str:
        logging.info("Missing static files generated, please update the following lines in the static-data section of your config file (written to stdout):")
        sys.stdout.write(update_str + '\n')
    else:
        logging.info("Missing static files generated.")

    return True
