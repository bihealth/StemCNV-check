import os
import sys
import tempfile
import ruamel.yaml as ruamel_yaml
from pathlib import Path

from deepdiff import DeepDiff
from pydantic.v1.utils import deep_update
from snakemake.api import SnakemakeApi
from snakemake.settings.types import ResourceSettings, ConfigSettings, DeploymentSettings, DAGSettings, OutputSettings, DeploymentMethod
from .check_input import check_config, check_sample_table
from stemcnv_check import STEM_CNV_CHECK, helpers
from loguru import logger as logging

import importlib.resources

def create_missing_staticdata(args):

    check_config(args, minimal_entries_only=True)
    config = helpers.load_config(args)
    check_sample_table(args)
    sample_df = helpers.read_sample_table(args.sample_table, args.column_remove_regex)

    ret = 0
    for array_name in config['array_definition']:
        # Check that at least one sample per array is defined
        # Otherwise give error, but continue with other arrays
        if not sample_df[sample_df['Array_Name'] == array_name].empty:
            ret += run_staticdata_workflow(args, array_name)
        else:
            logging.error(f'No samples defined for array "{array_name}"')
            ret += 1

    return ret


def run_staticdata_workflow(args, array_name):
    """
    Use the static_creation_global_files.smk workflow to generate missing static and array definition files.
    This includes download of (static) fasta and gtf files unless they are already defined in the config.
    For generation of array defnition files, the workflow runs the first steps of the main StemCNV-check 
    workflow to generate a vcf file if none are present. From this vcf file full information about GT
    (background) counts from the egt clusterfile can be extracted, that can be used to make the pfb file for PennCNV.
    Further array summary info is taken from UCSC genome info file and the probes locations from the vcf/PFB file.
    """
    logging.info(f'Starting to check for missing static data files for array "{array_name}" ...')
    # bpm & egt files are needed
    config = helpers.load_config(args)
    cache_path = helpers.get_cache_dir(args, config)

    genome_build = config['array_definition'][array_name]['genome_version']
    genome_build = 'hg38' if genome_build in ('hg38', 'GRCh38') else 'hg19'

    #short helper function for array files
    def array_file(filekey):
        return helpers.get_array_file(filekey, array_name, config, cache_path)

    # Collect missing global files
    missing_global_files = set()
    for global_file in ('fasta', 'gtf', 'genome_info', 'mehari_txdb', 'dosage_scores'):
        filename = helpers.get_global_file(global_file, genome_build, config['global_settings'], cache_path)
        if not os.path.isfile(filename):
            if global_file == 'mehari_txdb':
                os.makedirs(os.path.dirname(filename), exist_ok=True)
            missing_global_files.add(filename)

    if missing_global_files:
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
                        # nodes=args.jobs,
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
                        #apptainer_args=helpers.make_apptainer_args(config, cache_path, not_existing_ok=True),
                        conda_prefix=cache_path,
                        #apptainer_prefix=cache_path,
                    ),
                )
                .dag(
                    DAGSettings(targets=missing_global_files,
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
    static_file_keys = {'penncnv_pfb_file', 'penncnv_GCmodel_file', 'array_density_file', 'array_gaps_file'}
    # Update static_snake_config with existing and target file paths
    existing_files = {
        file for file in static_file_keys if file in config['array_definition'][array_name] and
        os.path.isfile(array_file(file))
    }
    missing_files = static_file_keys - existing_files
    static_snake_config.update({
        file: array_file(file) for file in existing_files
    })
    static_snake_config.update({
        file: array_file(file).format(
            genome=genome_build, array_name=array_name, cache=cache_path
        ) for file in missing_files
    })
    static_snake_config['generate_missing'] = list(missing_files)

    if not missing_files:
        logging.info('All static files are present')
        return 0

    # Check if vcf file _for the selected array_ is present, generate one if none are
    sample_data_df = helpers.read_sample_table(args.sample_table, args.column_remove_regex)
    sample_data_df = sample_data_df[sample_data_df['Array_Name'] == array_name]
    datapath = config['data_path']
    filter_settings = config['settings']['default_probe_filter_set']
    vcf_files = [
        os.path.join(datapath, f"{sample_id}", f"{sample_id}.annotated-SNP-data.{filter_settings}-filter.vcf.gz")
        for sample_id in sample_data_df['Sample_ID']
    ]
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
                            # nodes=args.jobs,
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
                            apptainer_args=helpers.make_apptainer_args(
                                config, cache_path, not_existing_ok=True, extra_bind_args=args.bind_points
                            ),
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
                    # nodes=args.jobs,
                ),
                config_settings=ConfigSettings(
                    config=dict(
                        static_snake_config,
                        **{
                            'TMPDIR': tmpdir,
                            'vcf_input_file': use_vcf,
                            'array_name': array_name,
                            'genomeInfo': helpers.get_global_file('genome_info', genome_build, config['global_settings'], cache_path),
                        }
                    )
                ),
                deployment_settings=DeploymentSettings(
                    deployment_method=DeploymentMethod.parse_choices_set(['conda', 'apptainer']),
                    conda_prefix=cache_path,
                    apptainer_prefix=cache_path,
                    apptainer_args=helpers.make_apptainer_args(
                        config, cache_path, tmpdir=tmpdir, not_existing_ok=True, extra_bind_args=args.bind_points
                    ),
                )
            )
            .dag(
                DAGSettings(targets={static_snake_config[file] for file in missing_files},
                            force_incomplete=True)
            )
            .execute_workflow()
        )

    # The missing files may have been specified with directly usable paths in the main config
    # If either format-strings or '__default-case__' was used the main config will differ from static_snake_config
    update_str = '\n'.join([
        f"    {entry}: {file}" for entry, file in static_snake_config.items() if
        entry in missing_files and
        entry in config['array_definition'][array_name] and
        file not in config['array_definition'][array_name][entry]
    ])
    if not update_str:
        logging.info(
            "Missing static files generated. The files have been written to the specified file paths and will not "
            "be written to a global array definition config."
        )
        return 0

    # Build full array config to write out: start from defined config
    # Then update all non-required static files, which may had format strings or '__default-case__' in the main config
    config_out = {'array_definition': {array_name: config['array_definition'][array_name]}}
    for file in static_file_keys:
        config_out['array_definition'][array_name][file] = static_snake_config[file]
    # In all cases the file paths should be absolute,
    # otherwise reusing that global config part from another path would break
    config_out['array_definition'][array_name] = {
        entry: os.path.abspath(value) if entry.endswith('_file') else value
        for entry, value in config_out['array_definition'][array_name].items()
    }

    # Determine where to write config updates
    # 1) if cache is in use, check if
    #    a) global config already exists, if not write to it
    #    b) check if the array definition already exists there, if not update global config
    #    c) if a (conflicting) array definition already exists, ask user if they want to overwrite
    # 2) if no cache is used, write to a local array specific config (with only that entry)
    #    Give an example for how to refer to this file in the main config file
    yaml = ruamel_yaml.YAML()
    base_config = {}
    array_name_sanitized = "".join(x if x.isalnum() or x in '_-.' else '_' for x in array_name)
    if cache_path:
        global_config_path = helpers.get_cache_array_definition(cache_path)
        if global_config_path:
            with open(global_config_path) as f:
                base_config = yaml.load(f)
            yaml_write_path = global_config_path
        else:
            yaml_write_path = os.path.join(cache_path, 'global_array_definitions.yaml')
    else:
        yaml_write_path = array_name_sanitized + '_config.yaml'

    if base_config and 'array_definition' in base_config and array_name in base_config['array_definition']:
        if base_config['array_definition'][array_name] == config_out['array_definition'][array_name]:
            logging.error("No changes to array definition detected, this should not happen after creating new files!")
            return 1
        # delta = DeepDiff(
        #     base_config['array_definition'][array_name],
        #     config_out['array_definition'][array_name],
        #     verbose_level=2
        # ).to_dict()
        logging.warning(
            f'Global config already contains array definition for "{array_name}".'
            # f' Differences:\n'
            # + '\n'.join([f"    {key}: {value}" for key, value in delta.items()])
        )
        answer = ''
        while not answer or answer[0] not in ('y', 'n'):
            answer = input(
                f'Do you want to overwrite the array definition for "{array_name}" in the gloabl config? (y/n): '
            ).lower()
        if answer == 'n':
            yaml_write_path = array_name_sanitized + '_config.yaml'
        else:
            config_out = deep_update(base_config, config_out)
    else:
        config_out = deep_update(base_config, config_out)

    info_str = (
        f'local file "{yaml_write_path}"' if (yaml_write_path == array_name_sanitized + '_config.yaml') else
        'global config.' #You either need to update the following entries in your main config ({args.config}) or remove '
        #f'the "{array_name}" array definition there, so that the global config can be used instead'
    )
    logging.info(
        f'Missing static files generated. Writing array definiton data entries to ' +
        info_str + ':\n' + update_str
    )
    with open(yaml_write_path, 'w') as f:
        yaml.dump(config_out, f)

    return 0
