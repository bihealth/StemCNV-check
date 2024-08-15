#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Wrapper around snakemake to kick-off pipeline
"""

import argparse
import os
import sys

from loguru import logger as logging

from . import __version__
from .app.check_input import check_config, check_sample_table
from .app.run_workflow import run_stemcnv_check_workflow
from .app.make_staticdata import create_missing_staticdata
from .app.setup_files import setup_control_files

def setup_argparse():

    parser = argparse.ArgumentParser(description="StemCNV-check: A pipeline to check the quality of SNP-array data for stem cell lines")

    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
    )

    parser.add_argument('action', nargs='?', default='run', choices=('run', 'setup-files', 'make-staticdata'),
                             help='Action to perform. Default: %(default)s')

    group_basic = parser.add_argument_group("General", "General pipeline arguments")

    group_basic.add_argument('--config', '-c', default='config.yaml', help="Filename of config file. Default: %(default)s")
    group_basic.add_argument('--sample-table', '-s', default='sample_table.tsv', help="Filename of sample table. Default: %(default)s")
    group_basic.add_argument('--directory', '-d', default=None,
                             help="Directory to run pipeline in. Default: current directory")

    # group_basic.add_argument('--non-interactive', '-y', action="store_true",
    #                          help="Skip all interactive questions of the wrapper")
    # group_basic.add_argument('--conda-frontend', default='mamba', choices=('mamba', 'conda'), help="Conda frontend to use. Default: %(default)s")
    group_basic.add_argument('--verbose', '-v', action='count', default=0,
                             help="More verbose output, maximum verbosity at -vv")

    group_setupfiles = parser.add_argument_group("setup-files", "Details for setup-files")
    group_setupfiles.add_argument('--config-details', default='minimal', choices=('minimal', 'medium', 'advanced', 'complete'), help="Level of detail for the config file. Default: %(default)s")
    group_setupfiles.add_argument('--overwrite', action='store_true', help="Allow overwriting of existing files")

    group_static = parser.add_argument_group("make-staticdata", "Details and file naming for make-staticdata")
    # group_static.add_argument('--genome', default='hg38', choices=('hg19', 'hg38'), help="Genome build to use (UCSC names). Default: %(default)s")
    # group_static.add_argument('--snp-array-name', default=None, help="A name or identifier string for the snp-array, can used in filesnames. No Default.")
    group_static.add_argument('--edit-config-inplace', action='store_true', help = "Edit the config file in place with updated static-data entries")
    group_static.add_argument('--penncnv-pfb-file', default='static-data/PennCNV-PFB_{genome}{array}.pfb',
                               help="Filename for generated PFB file. Default: %(default)s")
    group_static.add_argument('--penncnv-gcmodel-file', default='static-data/PennCNV-GCmodel_{genome}{array}.gcmodel',
                               help="Filename for generated GCmodel file. Default: %(default)s")
    group_static.add_argument('--array-density-file', default='static-data/density_{genome}{array}.bed',
                               help="Filename for generated bed file with probe density. Default: %(default)s")
    group_static.add_argument('--array-gaps-file', default='static-data/gaps_{genome}{array}.bed',
                              help="Filename for generated bed file with probe gaps. Default: %(default)s")
    group_static.add_argument('--genomeinfo-file', default='static-data/UCSC_{genome}_chromosome-info.tsv',
                               help="Filename for generated chromosome info file. Default: %(default)s")
    group_static.add_argument('--genome-gtf-file', default='static-data/gencode.{genome}.v45.gtf',
                               help="Filename for generated chromosome info file. Default: %(default)s")
    group_static.add_argument('--genome-fasta-file', default='static-data/{genome}.genome.fa.gz',
                               help="Filename for generated chromosome info file. Default: %(default)s")

    group_snake = parser.add_argument_group("Snakemake Settings", "Arguments for Snakemake (also affects make-staticdata)")

    group_snake.add_argument('--cache-path', default=None,
                             help="Override auto-selection of the cache path to a specific directory."
                             )
    group_snake.add_argument('--no-cache', action='store_true',
                             help="Do not use the a chache directory. The cache is used for workflow created metadata "
                             "(conda envs, singularity images, and VEP data). The default cache path is ~/.")
    
    group_snake.add_argument('--target', '-t', default='complete',
                             choices=('complete', 'report', 'combined-cnv-calls', 'PennCNV', 'CBS', 'SNP-data'),
                             help="Final target of the pipeline. Default: %(default)s")
    group_snake.add_argument('--cluster-profile', '-p', nargs='?', const='cubi-dev', help="Use snakemake profile for job submission to cluster. Default if used: %(const)s")
    group_snake.add_argument('-jobs', '-j', default=20, help="Number of oarallel job submissions in cluster mode. Default: %(default)s")
    group_snake.add_argument('--local-cores', '-n', default=4, help="Number of cores for local submission. Default: %(default)s")
    group_snake.add_argument('snake_options', nargs='*', #argparse.REMAINDER,
                             help="Options to pass to snakemake; separate from normal options with '--'")

    return parser


def main(argv=None):

    parser = setup_argparse()
    args = parser.parse_args(argv)

    if not args.verbose:
        sys.tracebacklimit = 0
    logging.remove(0)
    logging.add(sys.stderr,
                level=["WARNING", "INFO", "DEBUG"][min(args.verbose, 2)],
                backtrace=args.verbose > 0,
                diagnose=args.verbose > 1,
                )

    if args.action == 'run':
        check_sample_table(args.sample_table, args.config)
        check_config(args.config, args.sample_table)
        if args.directory is not None and not os.path.isdir(args.directory):
            os.makedirs(args.directory)
        ret = run_stemcnv_check_workflow(args)
    elif args.action == 'setup-files':
        ret = setup_control_files(args)
    elif args.action == 'make-staticdata':
        # args.snp_array_name = ('_' + args.snp_array_name) if args.snp_array_name and not args.snp_array_name[0] in ".-_" else ""
        # args.penncnv_pfb_file = args.penncnv_pfb_file.format(genome=args.genome, array=args.snp_array_name)
        # args.penncnv_gcmodel_file = args.penncnv_gcmodel_file.format(genome=args.genome, array=args.snp_array_name)
        # args.array_density_file = args.array_density_file.format(genome=args.genome, array=args.snp_array_name)
        # args.array_gaps_file = args.array_gaps_file.format(genome=args.genome, array=args.snp_array_name)
        # args.genomeinfo_file = args.genomeinfo_file.format(genome=args.genome)
        # args.genome_gtf_file = args.genome_gtf_file.format(genome=args.genome)
        # args.genome_fasta_file = args.genome_fasta_file.format(genome=args.genome)
        if args.directory is not None and not os.path.isdir(args.directory):
            os.makedirs(args.directory)
        ret = create_missing_staticdata(args)

    return ret


if __name__ == '__main__':
    sys.exit(main())

