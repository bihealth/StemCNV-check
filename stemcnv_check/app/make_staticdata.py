import os
import sys
import tempfile
import ruamel.yaml as ruamel_yaml
from pathlib import Path

from snakemake.api import SnakemakeApi
from snakemake.settings.types import ResourceSettings, ConfigSettings, DeploymentSettings, DAGSettings, OutputSettings, DeploymentMethod
from .check_input import check_config
from .. import STEM_CNV_CHECK
from ..helpers import config_extract, make_apptainer_args, read_sample_table, get_cache_dir, get_vep_cache_path
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

    # typ = 'safe' prevents round-trip writeout
    yaml = ruamel_yaml.YAML()
    with open(args.config) as f:
        config = yaml.load(f)

    config_def_file = importlib.resources.files(STEM_CNV_CHECK).joinpath('control_files', 'default_config.yaml')
    with config_def_file.open() as f:
        default_config = yaml.load(f)

    cache_path = get_cache_dir(args)
    use_cache = cache_path is not None

    vep_cache_path = get_vep_cache_path(
        config_extract(('settings', 'VEP_annotation', 'VEP_cache_path'), config, default_config),
        cache_path
    )
    get_vep_cache = (config_extract(('settings', 'VEP_annotation', 'mode'), config, default_config) == 'offline-cache' and 
                     config_extract(('settings', 'VEP_annotation', 'enabled'), config, default_config))

    array_name = config_extract(('array_name', ), config, default_config)
    genome_build = config_extract(('genome_version', ), config, default_config)
    genome_build = 'hg38' if genome_build in ('hg38', 'GRCh38') else 'hg19'
    vep_genome = 'GRCh38' if genome_build == 'hg38' else 'GRCh37'
    static_snake_config = {
        'TMPDIR': '',
        'genome': genome_build,
        'vcf_input_file': '',
        'density_windows': config_extract(('settings', 'array_attribute_summary', 'density.windows',), config, default_config),
        'min_gap_size': config_extract(('settings', 'array_attribute_summary', 'min.gap.size',), config, default_config),
        'vep_cache_path': vep_cache_path,
        'vep_cache': os.path.join(vep_cache_path, f'.{vep_genome}.done'),
        'vep_fasta_path': os.path.join(vep_cache_path, 'fasta'),
    }

    static_files = {
        #'genome_fasta_file', 
        'genome_gtf_file',
        'penncnv_pfb_file', 'penncnv_GCmodel_file',
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
    if get_vep_cache and not os.path.isfile(os.path.join(vep_cache_path, f'.{vep_genome}.done')):
        missing_files.add('vep_cache')

    genome_fasta = config_extract(('global_settings', f'{genome_build}_genome_fasta'), config, default_config)
    if genome_fasta == '__use-vep__':
        genome_fasta = os.path.join(vep_cache_path, 'fasta',
                                    'homo_sapiens', f'112_{vep_genome}',
                                    'Homo_sapiens.GRCh38.dna.toplevel.fa.gz' if vep_genome == 'GRCh38' else 
                                    'Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz'
                                    )
        get_fasta = not os.path.isfile(genome_fasta)
    else:
        get_fasta = False

    if not missing_files and not get_fasta:
        logging.info('All static files are present')
        return 0

    #TODO: will need annotated vcf soon, which also requires teh VEP cache to be present
    if get_fasta:
        logging.info(f'Genome fasta file for {genome_build} missing, will be downloaded with vep')
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
                        config=dict(static_snake_config, **{'TMPDIR': tmpdir, 'vep_fasta_path': vep_fasta_path})
                    ),
                    deployment_settings=DeploymentSettings(
                        deployment_method=DeploymentMethod.parse_choices_set({'conda', 'apptainer'}),
                        apptainer_args=make_apptainer_args(config, not_existing_ok=True),
                        conda_prefix=cache_path if use_cache else None,
                        apptainer_prefix=cache_path if use_cache else None,
                    ),
                )
                .dag(
                    DAGSettings(targets=[genome_fasta],
                                force_incomplete=True)
                )
                .execute_workflow()
            )
        # if not result:
        #     logging.error('Snakemake run to get fasta failed')
        #     sys.exit(1)
        # if fix_fasta_entry and not args.edit_config_inplace:
        #     logging.info(f"Please update the genome_fastq entry in the config and the restart this command.\n  genome_fasta_file: {args.genome_fasta_file}")
        #     sys.exit(0)
        # elif fix_fasta_entry:
        #     logging.info(f"Updating config file with new genome_fasta_file entry: {static_snake_config['genome_fasta_file']}")
        #     config['static_data']['genome_fasta_file'] = static_snake_config['genome_fasta_file']
        #     with open(args.config, 'w') as f:
        #         yaml.dump(config, f)

    # Check if vcf file is present, generate one if none are
    sample_data = read_sample_table(args.sample_table)
    datapath = config_extract(('data_path',), config, default_config)
    filter_settings = config_extract(('settings', 'default-filter-set'), config, default_config)
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
                                Path(args.config)
                            ],
                            config={
                                'sample_table': args.sample_table,
                                'basedir': args.directory,
                                'configfile': args.config,
                                'target': 'SNP-probe-data',
                                'use_vep_cache': vep_cache_path,
                                'global_settings': {f'{genome_build}_genome_fasta': genome_fasta}
                            }
                        ),
                        deployment_settings=DeploymentSettings(
                            deployment_method=DeploymentMethod.parse_choices_set(['conda', 'apptainer']),
                            conda_prefix=cache_path if use_cache else None,
                            apptainer_prefix=cache_path if use_cache else None,
                            apptainer_args=make_apptainer_args(config, not_existing_ok=True),
                        ),
                    )
                    .dag(
                        DAGSettings(targets=[use_vcf],
                                    force_incomplete=True)
                    )
                    .execute_workflow()
                )

        # if not result:
        #     logging.error('Snakemake run to get vcf failed')
        #     sys.exit(1)

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
                    conda_prefix=cache_path if use_cache else None,
                    apptainer_prefix=cache_path if use_cache else None,
                    apptainer_args=make_apptainer_args(config, tmpdir=tmpdir, not_existing_ok=True),
                )
            )
            .dag(
                DAGSettings(targets=[static_snake_config[file] for file in missing_files],
                            force_incomplete=True)
            )
            .execute_workflow()
        )

    update_str = '\n'.join([f"  {entry}: {file}" for entry, file in static_snake_config.items()
                            if entry in missing_files and entry != 'vep_cache' and
                            entry in config['static_data'] and file != config['static_data'][entry]])

    if args.edit_config_inplace and update_str:
        logging.info('Missing static files generated. Updating config file with new static data entries:\n' + update_str)
        for file in missing_files:
            if file == 'vep_cache':
                continue
            config['static_data'][file] = static_snake_config[file]
        with open(args.config, 'w') as f:
            yaml.dump(config, f)
    elif update_str:
        logging.info("Missing static files generated, please update the following lines in the static-data section of your config file (written to stdout):")
        sys.stdout.write(update_str + '\n')
    else:
        logging.info("Missing static files generated.")

    return True
