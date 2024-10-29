#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Wrapper around snakemake to kick-off pipeline
"""

import argparse
import datetime
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
    group_basic.add_argument('--sample-table', '-s', default=None,
                             help="Filename of sample table, can be tsv or xlsx format (1st sheet is read). "
                                  "Default: sample_table.tsv or sample_table.xlsx")
    group_basic.add_argument('--column-remove-regex', nargs='?', const=r'\s.*', type=str,
                             help="Regex to remove text from sample table column names (before looking for required columns)."
                                  " Not used by default, default if not regex given ' .*' (remove spaces and everything following a space)")
    # group_basic.add_argument('--required-col-indices', default=None, nargs=6,
    #                             help="Indexes of the required columns in the sample table. "
    #                                  "Required columns are: Sample_ID, Chip_Name, Chip_Pos, Array_Name, Sex, Reference_Sample. "
    #                                  "This will set the corresponding names to the specified columns")
    group_basic.add_argument('--directory', '-d', default=None,
                             help="Directory to run pipeline in. Default: current directory")
    group_basic.add_argument('--verbose', '-v', action='count', default=0,
                             help="More verbose output, maximum verbosity at -vv")

    group_setupfiles = parser.add_argument_group("setup-files", "Details for setup-files")
    group_setupfiles.add_argument('--config-details', default='minimal', choices=('minimal', 'medium', 'advanced', 'complete'), help="Level of detail for the config file. Default: %(default)s")
    group_setupfiles.add_argument('--sampletable-format', default='tsv', choices=('tsv', 'xlsx'), help="Format of the sample table. Default: %(default)s")
    group_setupfiles.add_argument('--overwrite', action='store_true', help="Allow overwriting of existing files")

    group_static = parser.add_argument_group("make-staticdata", "Details and file naming for make-staticdata")
    group_static.add_argument('--no-edit-inplace', action='store_true', help = "Do not edit the config file in place with updated static-data entries")

    group_snake = parser.add_argument_group("Snakemake Settings", "Arguments for Snakemake (also affects make-staticdata)")
    group_snake.add_argument('--cache-path', default=None,
                             help="Override auto-selection of the cache path to a specific directory. The default cache path is defined in the conifg file."
                             )
    group_snake.add_argument('--no-cache', action='store_true',
                             help="Do not use a chache directory. The cache is used for workflow created metadata "
                             "(conda envs, singularity images, and VEP data). The default cache path is defined in the conifg file.")

    group_snake.add_argument('--target', '-t', default='complete',
                             choices=('complete', 'report', 'collate-summary', 'summary-tables', 'collate-cnv-calls',
                                      'combined-cnv-calls', 'PennCNV', 'CBS', 'SNP-data'),
                             help="Final target of the pipeline. Default: %(default)s")
    group_snake.add_argument('--collate-date', nargs='?', const=datetime.date.today().strftime("%Y-%m-%d"),
                             default=None, help="Add a date to the collate output files. Default without argument: today's date")

    group_snake.add_argument('--cluster-profile', '-p', help="Use snakemake profile for job submission to cluster. Default if used: %(const)s")
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

    if args.sample_table is None and args.action != 'setup-files':
        if os.path.isfile('sample_table.tsv'):
            args.sample_table = 'sample_table.tsv'
        elif os.path.isfile('sample_table.xlsx'):
            args.sample_table = 'sample_table.xlsx'
        else:
            logging.error("No default sample table found (sample_table.tsv or sample_table.xlsx). "
                          "Please create a sample table (i.e. stemcnv-check setup-files) or specify one with --sample-table")
            raise FileNotFoundError("No sample table found")

    if args.action == 'run':
        check_sample_table(args.sample_table, args.config, args.column_remove_regex)
        check_config(args.config, args.sample_table, args.column_remove_regex)
        if args.directory is not None and not os.path.isdir(args.directory):
            os.makedirs(args.directory)
        ret = run_stemcnv_check_workflow(args)
    elif args.action == 'setup-files':
        if not args.sample_table:
            args.sample_table = 'sample_table.tsv' if args.sampletable_format == 'tsv' else 'sample_table.xlsx'
        ret = setup_control_files(args)
    elif args.action == 'make-staticdata':
        if args.directory is not None and not os.path.isdir(args.directory):
            os.makedirs(args.directory)
        ret = create_missing_staticdata(args)

    return ret


if __name__ == '__main__':
    sys.exit(main())

