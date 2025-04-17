#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Wrapper around snakemake to kick-off pipeline
"""

import argparse
import datetime
import multiprocessing
import os
import platform
import psutil
import sys

from loguru import logger as logging

from . import __version__
from .app.check_input import check_config, check_sample_table
from .app.run_workflow import run_stemcnv_check_workflow
from .app.make_staticdata import create_missing_staticdata
from .app.setup_files import setup_control_files

# from: https://github.com/scivision/detect-windows-subsystem-for-linux/blob/main/is_wsl.py
def is_wsl(v: str = platform.uname().release) -> int:
    """
    detects if Python is running in WSL
    """

    if v.endswith("-Microsoft"):
        return 1
    elif v.endswith("microsoft-standard-WSL2"):
        return 2

    return 0

def get_default_memory_limit():

    if not is_wsl():
        return None
    # available memory in MB
    mem_mb = psutil.virtual_memory().total / (1024 * 1024)
    return int(mem_mb)

def setup_argparse():
    # Check if running on WSL, if yes we need to set memory limit
    
    # General options, available through inheritance for all subparsers
    general_options = argparse.ArgumentParser(description="General StemCNV-check options", add_help=False)
    general_group = general_options.add_argument_group("General", "General pipeline arguments")
    general_group.add_argument(
        '-v', '--verbose', action='count', default=0,
        help="More verbose output and additional logging, maximum verbosity at -vv"
    )
    general_group.add_argument(
        '-c', '--config', default='config.yaml',
        help="Filename of config file. Default: %(default)s\n"
             "Note: if a global array definition exists in the cache path, it will also be used by default"
    )
    general_group.add_argument(
         '-s', '--sample-table', default=None,
        help="Filename of sample table, can be tsv or xlsx format (1st sheet is read). "
             "Default: sample_table.tsv or sample_table.xlsx"
    )

    # Base parser
    parser = argparse.ArgumentParser(
        description="StemCNV-check: A pipeline to check the quality of SNP-array data for stem cell lines",
    )

    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
    )
    subparsers = parser.add_subparsers(
        dest='action', required=True,
    )

    group_setupfiles = subparsers.add_parser(
        "setup-files",
        help="Create a sample table and config file for the pipeline",
        parents=[general_group],
    )
    group_setupfiles.add_argument(
        '--config-details', default='minimal',
        choices=('minimal', 'medium', 'advanced', 'complete'),
        help="Level of detail for the config file. Default: %(default)s"
    )
    group_setupfiles.add_argument(
        '--sampletable-format', default='tsv',
        choices=('tsv', 'xlsx'),
        help="Format of the sample table. Default: %(default)s"
    )
    group_setupfiles.add_argument(
        '--overwrite', action='store_true',
        help="Allow overwriting of existing files"
    )

    # Snakemake options, used and inherited by both run & make-staticdata 
    group_snake = argparse.ArgumentParser(description="StemCNV-check snakemake options", add_help=False)

    group_snake.add_argument(
        '-d', '--directory', default=None,
        help="Directory to run pipeline in. Default: current directory"
    )
    group_snake.add_argument(
        '-n', '--local-cores', default=multiprocessing.cpu_count(),
        help="Number of cores for local submission. Default is to use all available (%(default)s)"
    )
    group_snake.add_argument(
        '--memory-mb', default=get_default_memory_limit(),
        help="Maximum memory to use, in Mb. Default is None on Linux, on WSL the detected available memory: %(default)s"
    )
    group_snake.add_argument(
        '--cache-path', default=None,
        help="Override auto-selection of the cache path to a specific directory."
             " The default cache path is defined in the conifg file."
    )
    group_snake.add_argument(
        '--no-cache', action='store_true',
        help="Do not use a cache directory for workflow created metadata. (cache includes: global array definition "
             "config, conda envs, singularity images, and reference data). "
             "The default cache path is defined in the conifg file."
    )
    group_snake.add_argument(
        '--column-remove-regex', nargs='?', const=r'\s.*', type=str,
        help="Regex to remove text from sample table column names (before looking for required columns)."
             " Not used by default, default if no regex given is ' .*' (remove spaces and everything following a space)"
    )
    group_snake.add_argument(
        '--bind-points',
        help="Additional bind points for apptainer containers, intended for expert users. Use i.e. '/path' to make "
             "it available in apptainer, useful in case local directory "
             "contains symlinks that won't resolve in the container."
    )

    group_staticdata = subparsers.add_parser(
        "make-staticdata", help="Create StemCNV-check static data files using snakemake",
        parents=[general_group, group_snake]
    )


    group_run = subparsers.add_parser(
        "run", #dest="run",
        help="Run the StemCNV-check snakemake workflow",
        parents=[general_group, group_snake]
    )

    group_run.add_argument(
        '--target', '-t', default='complete',
        choices=(
            'complete', 'report', 'collate-summary', 'summary-tables', 'collate-cnv-calls',
            'combined-cnv-calls', 'PennCNV', 'CBS', 'SNP-data', 'gtc-data'
        ),
        help="Final target of the pipeline. Default: %(default)s"
    )
    group_run.add_argument(
        '--collate-date', default=None, nargs='?',
        const=datetime.date.today().strftime("%Y-%m-%d"),
        help="Add a date to the collate output files. Default without argument: today's date"
    )
    group_run.add_argument(
        '--snakemake-help', action='store_true',
        help="Show snakemake help message, all snakemake options must be passed after '--'"
    )
    group_run.add_argument('snake_options', nargs='*', #argparse.REMAINDER,
                             help="Options to pass to snakemake; separate from normal options with '--'")

    return parser


def main(argv=None):

    parser = setup_argparse()
    args = parser.parse_args(argv)

    if args.action == 'run' and args.snakemake_help:
        from snakemake.cli import main as snakemake_main
        snakemake_main(['--help'])
        return 0

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
        check_sample_table(args)
        check_config(args)
        if args.directory is not None and not os.path.isdir(args.directory):
            os.makedirs(args.directory)
        ret = run_stemcnv_check_workflow(args, is_wsl())
    elif args.action == 'setup-files':
        if not args.sample_table:
            args.sample_table = 'sample_table.tsv' if args.sampletable_format == 'tsv' else 'sample_table.xlsx'
        ret = setup_control_files(args)
    elif args.action == 'make-staticdata':
        if args.directory is not None and not os.path.isdir(args.directory):
            os.makedirs(args.directory)
        ret = create_missing_staticdata(args, is_wsl())

    return ret


if __name__ == '__main__':
    sys.exit(main())

