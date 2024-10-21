import os
import sys
import tempfile
import ruamel.yaml as ruamel_yaml
from pathlib import Path

from snakemake.api import SnakemakeApi
from snakemake.settings.types import ResourceSettings, ConfigSettings, DeploymentSettings, DAGSettings, OutputSettings, DeploymentMethod
from .check_input import check_config
from .. import STEM_CNV_CHECK
from ..helpers import make_apptainer_args, read_sample_table, get_cache_dir, load_config, get_global_file
from loguru import logger as logging

import importlib.resources

def create_missing_staticdata(args):

    check_config(args.config, args.sample_table, args.column_remove_regex, required_only=True)
    config = load_config(args.config)

    for array_name in config['array_definition']:
        if array_name == '_default_':
            continue
        run_staticdata_workflow(args, array_name)


def run_staticdata_workflow(args, array_name):
    """
    Use the static_creation_global_files.smk workflow to generate missing static and array definition files.
    This includes download of (static) fasta and gtf files unless they are already defined in the config.
    For generation of array defnition files, the workflow runs the first steps of the main StemCNV-check 
    workflow to generate a vcf file if none are present. From this vcf file full information about GT
    (background) counts from the egt clusterfile can be extracted, that can be used to make the pfb file for PennCNV.
    Further array summary info is take from UCSC genome info file and the probes locations from the vcf/PFB file.
    """
    logging.info(f'Starting to check for missing static data files for array "{array_name}" ...')
    # bpm & egt files are needed
    config = load_config(args.config)
    cache_path = get_cache_dir(args, config)

    genome_build = config['array_definition'][array_name]['genome_version']
    genome_build = 'hg38' if genome_build in ('hg38', 'GRCh38') else 'hg19'

    # Collect missing global files
    global_files = set()
    for global_file in ('fasta', 'gtf', 'genome_info', 'mehari_txdb'):
        filename = get_global_file(global_file, genome_build, config['global_settings'], cache_path)
        if not os.path.isfile(filename):
            if global_file == 'mehari_txdb':
                os.makedirs(os.path.dirname(filename), exist_ok=True)
            global_files.add(filename)

    if global_files:
        logging.info(f'Setting up genome wide static files {genome_build}')
        # This is now a separated workflow, with doesn't need apptainer (and makes the basic files always bound to apptainer)
        with (SnakemakeApi(output_settings=OutputSettings(printshellcmds=True)) as api,
              tempfile.TemporaryDirectory() as tmpdir):
            result = (
                api.workflow(
                    snakefile=importlib.resources.files(STEM_CNV_CHECK).joinpath('rules', 'static_creation_global_files.smk'),
                    workdir=args.directory,
                    resource_settings=ResourceSettings(
                        cores=args.local_cores,
                        local_cores=args.local_cores,
                        nodes=args.jobs,
                    ),
                    config_settings=ConfigSettings(
                        config={
                            'global_settings': config['global_settings'],
                            'genome': genome_build,
                            'cache_path': cache_path,
                            'TMPDIR': tmpdir
                        }
                    ),
                    deployment_settings=DeploymentSettings(
                        deployment_method=DeploymentMethod.parse_choices_set(['conda']),
                        #apptainer_args=make_apptainer_args(config, cache_path, not_existing_ok=True),
                        conda_prefix=cache_path,
                        #apptainer_prefix=cache_path,
                    ),
                )
                .dag(
                    DAGSettings(targets=global_files,
                                force_incomplete=True)
                )
                .execute_workflow()
            )

    static_snake_config = {
        'TMPDIR': '',
        'genome': genome_build,
        'vcf_input_file': '',
        'density_windows': config['settings']['array_attribute_summary']['density.windows'],
        'min_gap_size': config['settings']['array_attribute_summary']['min.gap.size'],
    }
    # Required auto-generatbale array definition files
    static_files = {'penncnv_pfb_file', 'penncnv_GCmodel_file', 'array_density_file', 'array_gaps_file'}
    # Update static_snake_config with existing and target file paths
    existing_files = {file for file in static_files
                      if file in config['array_definition'][array_name] and os.path.isfile(config['array_definition'][array_name][file])}
    missing_files = static_files - existing_files
    static_snake_config.update({
        file: config['array_definition'][array_name][file] for file in existing_files
    })
    static_snake_config.update({
        file: config['array_definition'][array_name][file].format(genome=genome_build, array_name=array_name) for file in missing_files
    })

    if not missing_files:
        logging.info('All static files are present')
        return 0

    # Check if vcf file is present, generate one if none are
    sample_data = read_sample_table(args.sample_table, args.column_remove_regex)
    datapath = config['data_path']
    filter_settings = config['settings']['default-filter-set']
    #FIXME (future): check if annotation is enabled/disabled
    vcf_files = [os.path.join(datapath, f"{sample_id}", f"{sample_id}.annotated-SNP-data.{filter_settings}-filter.vcf.gz") for
                 sample_id, _, _, _, _, _ in sample_data]
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
                                'column_remove_regex': args.column_remove_regex,
                                'basedir': args.directory,
                                # 'configfile': args.config,
                                # 'target': 'SNP-probe-data',
                                'cache_path': str(cache_path),
                                'verbose_level': args.verbose,
                            }
                        ),
                        deployment_settings=DeploymentSettings(
                            deployment_method=DeploymentMethod.parse_choices_set(['conda', 'apptainer']),
                            conda_prefix=cache_path,
                            apptainer_prefix=cache_path,
                            apptainer_args=make_apptainer_args(config, cache_path, not_existing_ok=True),
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
                snakefile=importlib.resources.files(STEM_CNV_CHECK).joinpath('rules', 'static_creation_array_definition.smk'),
                workdir=args.directory,
                resource_settings=ResourceSettings(
                    cores=args.local_cores,
                    local_cores=args.local_cores,
                    nodes=args.jobs,
                ),
                config_settings=ConfigSettings(
                    config=dict(
                        static_snake_config,
                        **{
                            'TMPDIR': tmpdir,
                            'vcf_input_file': use_vcf,
                            'genomeInfo': get_global_file('genome_info', genome_build, config['global_settings'], cache_path),
                        }
                    )
                ),
                deployment_settings=DeploymentSettings(
                    deployment_method=DeploymentMethod.parse_choices_set(['conda', 'apptainer']),
                    conda_prefix=cache_path,
                    apptainer_prefix=cache_path,
                    apptainer_args=make_apptainer_args(config, cache_path, tmpdir=tmpdir, not_existing_ok=True),
                )
            )
            .dag(
                DAGSettings(targets={static_snake_config[file] for file in missing_files},
                            force_incomplete=True)
            )
            .execute_workflow()
        )

    update_str = '\n'.join([f"    {entry}: {file}" for entry, file in static_snake_config.items()
                            if entry in missing_files and
                            entry in config['array_definition'][array_name] and file != config['array_definition'][array_name][entry]])

    if not args.no_edit_inplace and update_str:
        logging.info('Missing static files generated. Updating config file with array definiton data entries:\n' + update_str)
        user_config = load_config(args.config, False)
        for file in missing_files:
            user_config['array_definition'][array_name][file] = static_snake_config[file]
        yaml = ruamel_yaml.YAML()
        with open(args.config, 'w') as f:
            yaml.dump(user_config, f)
    elif update_str:
        logging.info(f"Missing static files generated, please use the following lines to update the array_definition "
                     f"for '{array_name}' in your config file (written to stdout):")
        sys.stdout.write(update_str + '\n')
    else:
        logging.info("Missing static files generated.")

    return True
