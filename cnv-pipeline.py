#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Wrapper around snakemake to kick-off pipeline
"""

import argparse
#import datetime
#import logging
import os
import subprocess
import sys

# import ruamel.yaml as ruamel_yaml
# from snakemake import RERUN_TRIGGERS
from snakemake import main as snakemake_main
from scripts.py_helpers import read_sample_table


SNAKEDIR = os.path.dirname(os.path.realpath(__file__))

### Sanity checks ###

#TODO:
# sanity check of sample table (fail early)
# - all reference samples exist
# - gender is allowed values
# - FUTURE: sample_id / chip_No/Pos conversion/check
def check_sample_table(filename):
	pass

def check_config(filename):
	pass

def check_installation():
	pass


### Actions ###

def copy_setup_files(args):
	#TODO: copy into PWD or into directroy?
	subprocess.call(['cp', f"{SNAKEDIR}/config.yaml", args.config])
	subprocess.call(['cp', f"{SNAKEDIR}/sample_table.txt", args.sample_table])
	print('Created setup files: ...')

def make_penncnv_files(args):
	
	# Check if any vcf file is present
	sample_data = read_sample_table(args.sample_table)
	
	vcf_files = [os.path.join(args.directory, "data", f"{sample_id}", f"{sample_id}.unprocessed.vcf") for _, sample_id, _, _ in sample_data.values()]
	vcf_present = [vcf for vcf in vcf_files if os.path.exists(vcf)]
	
	if vcf_present:
		use_vcf = vcf_present[0]
	else:
		use_vcf = vcf_files[0]
		print('Running snakemake to get:', use_vcf)
		args.snake_options += [
			use_vcf
			#'--until', use_vcf,
			#'--allowed-rules', 'run_gtc2vcf_vcf', 
			#'run_gencall', 'relink_gencall', 'all'
		]
		ret = run_snakemake(args)
		
		if ret:
			raise Exception('Snakemake run to get vcf failed')
	
	#TODO logger calls
	print('Making PFB file')
	argv = ['Rscript', f"{SNAKEDIR}/scripts/make_PFB_from_vcf.R", use_vcf, args.pfb_out]
	print(' '.join(argv))
	ret = subprocess.call(argv)
	
	if ret:
		raise Exception('Generation of PFB from vcf file failed')
	else:
		print('Created PFB file: {}'.format(args.pfb_out))
	
	print('Making GC model file')
	gc_base = 'hg38.gc5Base.txt' if args.genome == 'GRCh38' else 'hg19.gc5Base.txt'
	penncnv_gcfile = os.path.join(os.getenv('CONDA_PREFIX'), 'pipeline/PennCNV-1.0.5/gc_file', gc_base)
	
	if not os.path.exists(penncnv_gcfile) and os.path.exists(penncnv_gcfile + '.gz'):
		subprocess.call(['gunzip', penncnv_gcfile + '.gz'])
	elif not os.path.exists(penncnv_gcfile):
		raise Exception('Could not find PennCNV GC files - did you activate the cnv-pipeline conda env?')
		
	argv = [os.path.join(os.getenv('CONDA_PREFIX'), 'pipeline/PennCNV-1.0.5/cal_gc_snp.pl'),
					penncnv_gcfile, args.pfb_out, '-out', args.gc_out]
	#print(' '.join(argv))
	ret = subprocess.call(argv)				
	
	if ret:
		raise Exception('Generation of GCmodel file with PennCNV')
	else:
		print('Created GCmodel file: {}'.format(args.gc_out))

	return 0


def run_snakemake(args):
	
	argv = [
		"-s", os.path.join(SNAKEDIR, "cnv-snakefile.py"),
		"-p", "-r",
		#"--keep-incomplete",
		"--rerun-incomplete"
	]
	
	argv += [
		'-d', args.directory,
		'--configfile', os.path.join(SNAKEDIR, 'default_config.yaml'), args.config,
		'--config', 'sample_table={}'.format(args.sample_table),
			'snakedir={}'.format(SNAKEDIR), 
			'basedir={}'.format(args.directory), 
			'configfile={}'.format(args.config),
	]
	
	if args.cluster:
		argv += [
			"--profile", "cubi-dev",
			"-j", "20"
		]
	else:
		argv += [
			'--cores', str(args.local_cores)
		]
	
	if args.snake_options:
		argv += args.snake_options
		
	#print(argv)
	
	return snakemake_main(argv)


### Setup ###

def setup_argparse():
	parser = argparse.ArgumentParser(description="run CNV pipeline and helpers")

	group_basic = parser.add_argument_group("General", "Genera pipeline arguments")

	group_basic.add_argument('--action', '-a', default= 'run', choices = ('run', 'make-penncnv-files'), help = 'Action to perform. Default: %(default)s')
	group_basic.add_argument('--config', default='config.yaml', help="Filename of config file. Default: %(default)s")
	group_basic.add_argument('--sample-table', '-s', default='sample_table.txt', help="Filename of sample table. Default: %(default)s")
	
	group_penncnv = parser.add_argument_group("make-penncnv-files", "Specific arguments for make-penncnv-files")
	group_penncnv.add_argument('--genome', default = 'GRCh38', choices = ('GRCh37', 'GRCh38'), 
														 help="Genome build to make the GC model for (uses the files shipped with PennCNV). Default: %(default)s")
	group_penncnv.add_argument('--pfb-out', default = 'static-data/PennCNV-PFB_from_clusterfile-stats.pfb',
														 help="Filename for generated PFB file. Default: %(default)s")
	group_penncnv.add_argument('--gc-out', default = 'static-data/PennCNV-GCmodel-GRCh38.gcmodel',
														 help="Filename for generated GCmodel file. Default: %(default)s")
														 
	group_snake = parser.add_argument_group("Snakemake Settings", "Arguments for Snakemake")
	
	group_snake.add_argument('--cluster', '-c', action='store_true', help="Use slurm submission to run on cluster")
	group_snake.add_argument('--local-cores', '-n', default=4, help="Number of cores for local submission. Default: %(default)s")
	group_snake.add_argument('--directory', '-d', default=os.getcwd(), help="Directory to run pipeline in. Default: $CWD")
	group_snake.add_argument('snake_options', nargs='*', #argparse.REMAINDER, 
														help="Options to pass to snakemake. Use seprate with '--'")
	
	return parser
	


if __name__ == '__main__':
	
	args = setup_argparse().parse_args()
	
	#TODO
	# check that conda & PennCNV are set up ?!
	
	if args.action == 'run':
		check_sample_table(args.sample_table)
		ret = run_snakemake(args)
	elif args.action == 'setup':
		ret = run_config_setup(args)
	elif args.action == 'make-penncnv-files':
		ret = make_penncnv_files(args)

	sys.exit(ret)

