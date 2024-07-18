import os
import sys
import tempfile
import ruamel.yaml as ruamel_yaml

from snakemake import snakemake
from .check_input import check_config
from .. import STEM_CNV_CHECK
from ..helpers import config_extract, make_singularity_args, read_sample_table
from loguru import logger as logging

import importlib.resources

def create_missing_staticdata(args):
    # Check if any vcf file is present, generate one if none are
    # This is done because the vcf files generated with gtc2vcf contain the full information
    # about GT (background) counts from the egt clusterfile, that can be used to make the pfb file for PennCNV
    logging.info('Starting to check for missing static data files ...')
    # bpm & egt files are needed
    check_config(args, required_only=True)

    # typ = 'safe' prevents round-trip writeout
    yaml = ruamel_yaml.YAML()
    with open(args.config) as f:
        config = yaml.load(f)

    default_config = importlib.resources.files(STEM_CNV_CHECK).joinpath('control_files', 'default_config.yaml')
    with default_config.open() as f:
        DEF_CONFIG = yaml.load(f)


    static_snake_config = {
        'snakedir': importlib.resources.files(STEM_CNV_CHECK),
        'use_singularity': not args.no_singularity,
        'TMPDIR': '',
        'genome': args.genome,
        'vcf_input_file': '',
        'density_windows': config_extract(('settings', 'array_attribute_summary', 'density.windows',), config, DEF_CONFIG),
        'min_gap_size': config_extract(('settings', 'array_attribute_summary', 'min.gap.size',), config, DEF_CONFIG),
        }

    static_files = ('genome_fasta_file', 'genome_gtf_file', 'penncnv_pfb_file', 'penncnv_GCmodel_file',
                 'genomeInfo_file', 'array_density_file', 'array_gaps_file')

    # fasta file needs to be available before the other files
    get_fasta = False
    fix_fasta_entry = False
    for file in static_files:
        if file not in config['static_data'] or not config['static_data'][file]:
            static_snake_config.update({file: getattr(args, file.lower())})
            if file == 'genome_fasta_file':
                get_fasta = True
                fix_fasta_entry = True
        elif not os.path.isfile(config['static_data'][file]):
            static_snake_config.update({file: getattr(args, file.lower())})
            if file == 'genome_fasta_file':
                get_fasta = True
                fix_fasta_entry = config['static_data'][file] != args.genome_fasta_file

    if get_fasta:
        logging.info('Genome fasta file not found in config or missing, will be downloaded from GenCode')
        with tempfile.TemporaryDirectory() as tmpdir:
            ret = snakemake(
                importlib.resources.files(STEM_CNV_CHECK).joinpath('rules', 'staticdata_creation.smk'),
                #os.path.join(codebase_dir, "staticdata_creation.smk"),
                local_cores=args.local_cores,
                cores=args.local_cores,
                workdir=args.directory,
                use_singularity=not args.no_singularity,
                singularity_args='' if args.no_singularity else make_singularity_args(config, not_existing_ok=True),
                use_conda=True,
                conda_frontend=args.conda_frontend,
                printshellcmds=True,
                force_incomplete=True,
                config=dict(static_snake_config, **{'TMPDIR': tmpdir, "static-data":{"genome_fasta_file": args.genome_fasta_file}}),
                targets=[args.genome_fasta_file]
            )
        if not ret:
            logging.error('Snakemake run to get fasta failed')
            sys.exit(1)
        if fix_fasta_entry and not args.edit_config_inplace:
            logging.info(f"Please update the genome_fastq entry in the config and the restart this command.\n  genome_fasta_file: {args.genome_fasta_file}")
            sys.exit(0)
        elif fix_fasta_entry:
            logging.info(f"Updating config file with new genome_fasta_file entry: {args.genome_fasta_file}")
            config['static_data']['genome_fasta_file'] = args.genome_fasta_file
            with open(args.config, 'w') as f:
                yaml.dump(config, f)


    # Check if vcf file is present, generate one if none are
    sample_data = read_sample_table(args.sample_table)
    datapath = config_extract(('data_path',), config, DEF_CONFIG)
    vcf_files = [os.path.join(args.directory, datapath, f"{sample_id}", f"{sample_id}.unprocessed.vcf") for
                 sample_id, _, _, _, _ in sample_data]
    vcf_present = [vcf for vcf in vcf_files if os.path.exists(vcf)]

    if vcf_present:
        use_vcf = vcf_present[0]
    else:
        use_vcf = vcf_files[0]
        logging.info(f'Running first steps of StemCNV-check to get a vcf file: {use_vcf}')
        ret = snakemake(
            importlib.resources.files(STEM_CNV_CHECK).joinpath('rules', 'StemCNV-check.smk'),
            local_cores=args.local_cores,
            cores=args.local_cores,
            workdir=args.directory,
            configfiles=[importlib.resources.files(STEM_CNV_CHECK).joinpath('control_files', 'default_config.yaml'), args.config],
            printshellcmds=True,
            force_incomplete=True,
            use_conda=True,
            use_singularity=not args.no_singularity,
            singularity_args='' if args.no_singularity else make_singularity_args(config, not_existing_ok=True),
            config={'sample_table': args.sample_table,
                    'snakedir': importlib.resources.files(STEM_CNV_CHECK),
                    'basedir': args.directory,
                    'configfile': args.config,
                    'use_singularity': not args.no_singularity
                    },
            targets=[use_vcf]
        )

        if not ret:
            logging.error('Snakemake run to get vcf failed')
            sys.exit(1)

    # Run extra snakemake to check which files are missing & create them accordingly
    logging.info(f'Running staticdata creation workflow')
    with tempfile.TemporaryDirectory() as tmpdir:
        ret = snakemake(
            importlib.resources.files(STEM_CNV_CHECK).joinpath('rules', 'staticdata_creation.smk'),
            local_cores=args.local_cores,
            cores=args.local_cores,
            workdir=args.directory,
            use_singularity=not args.no_singularity,
            singularity_args='' if args.no_singularity else make_singularity_args(config, tmpdir, True),
            use_conda=True,
            conda_frontend=args.conda_frontend,
            printshellcmds=True,
            force_incomplete=True,
            config=dict(static_snake_config, **{'TMPDIR': tmpdir, 'vcf_input_file': use_vcf}),
        )

    update_str = '\n'.join([f"  {entry}: {file}" for entry, file in static_snake_config.items() if entry in static_files])

    if ret and args.edit_config_inplace:
        logging.info('Updating config file with new static data entries:\n' + update_str)
        for file in static_files:
            if file in static_snake_config.keys():
                config['static_data'][file] = static_snake_config[file]
        with open(args.config, 'w') as f:
            yaml.dump(config, f)
    elif ret:
        logging.info("Missing static files generated, please update the following lines in the static-data section of your config file (written to stdout):")
        sys.stdout.write(update_str)

    return ret
