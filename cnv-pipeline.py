#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Wrapper around snakemake to kick-off pipeline
"""

import argparse
# import datetime
# import logging
import os
import re
import subprocess
import sys
import yaml

# import ruamel.yaml as ruamel_yaml
# from snakemake import RERUN_TRIGGERS
from snakemake import main as snakemake_main
from scripts.py_helpers import read_sample_table
from scripts.py_exceptions import *

SNAKEDIR = os.path.dirname(os.path.realpath(__file__))

### Sanity checks ###
def check_sample_table(args):
	sample_data = read_sample_table(args.sample_table)

	samples = {sid: sex for _, _, sid, sex, _ in sample_data.values()}
	ref_samples = {rid: sex for _, _, _, sex, rid in sample_data.values() if rid}

	#Check sex values
	if not all(s in ('m', 'f') for s in map(lambda x: x[0].lower(), samples.values())):
		raise SampleConstraintError("Not all values of the 'Sex' column in the samplesheet can be coerced to 'm' or 'f'")
	#Check that all reference samples exist
	missing_refs = [ref for ref in ref_samples.keys() if ref not in samples.keys()]
	if missing_refs:
		raise SampletableReferenceError("These 'Reference_Sample's do not also exist in the 'Sample_ID' column of the samplesheet: " + ', '.join(missing_refs))
	# Give warning if sex of reference and sample don't match
	sex_mismatch = [f"{s} ({sex})" for _, _, s, sex, ref in sample_data.values() if ref and sex[0].lower() != samples[ref][0].lower()]
	if sex_mismatch:
		sys.stderr.write("Warning: the following samples have a different sex annotation than their Reference_Sample: " + ', '.join(sex_mismatch))
	# Check that Chip_Name & Chip_Pos match the sentrix wildcard regex
	with open(args.config) as f:
		config = yaml.safe_load(f)
	with open(os.path.join(SNAKEDIR, 'default_config.yaml')) as f:
		def_config = yaml.safe_load(f)

	for constraint, val in (('sample_id', 's'), ('sentrix_name', 'n'), ('sentrix_pos', 'p')):
		pattern = config['wildcard_constraints'][constraint] if 'wildcard_constraints' in config and constraint in config['wildcard_constraints'] else def_config['wildcard_constraints'][constraint]
		mismatch = ["{} ({})".format(sn, eval(val)) for sn, (n, p, s, _, _) in sample_data.items() if not re.match('^' + pattern + '$', eval(val))]
		if mismatch:
			raise SampleConstraintError(f"The '{constraint}' values for these samples do not fit the expected constraints: " + ', '.join(mismatch))


# Assume that the default_config hasn't been altered & is correct
def check_config(args):

	with open(args.config) as f:
		config = yaml.safe_load(f)

	#TODO: need helper function to ignore values omitted from user config

	if config['settings']['make_cnv_vcf']['mode'] not in c('combined-calls', 'split-tools'):
		raise ConfigValueError(
			'Value not allowed for settings$make.cnv.vcf$mode: "{}"'.format(config['settings']['make.cnv.vcf']['mode']))

	#TODO:
	# - check that required values are filled in
	# - check that value types are correct / match config (incl list len?)
	# - maybe: us the same apprach that seasnap has for this?


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
	with open(args.config) as f:
		config = yaml.safe_load(f)
	with open(os.path.join(SNAKEDIR, 'default_config.yaml')) as f:
		def_config = yaml.safe_load(f)
	datapath = config['data_path'] if 'data_path' in config else def_config['data_path']
	
	vcf_files = [os.path.join(args.directory, datapath, f"{sample_id}", f"{sample_id}.unprocessed.vcf") for _, _, sample_id, _, _ in sample_data.values()]
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
		#TODO: why does the code terminate here (cluster mode only?)?
		
		if ret:
			raise Exception('Snakemake run to get vcf failed')
	
	#TODO logger calls
	print('Making PFB file')
	outdir = os.path.dirname(args.pfb_out)
	if not os.path.exists(outdir):
		os.makedirs(outdir)

	argv = ['Rscript', f"{SNAKEDIR}/scripts/make_PFB_from_vcf.R", use_vcf, args.pfb_out]
	print(' '.join(argv))
	ret = subprocess.call(argv)
	
	if ret:
		raise Exception('Generation of PFB from vcf file failed')
	else:
		print('Created PFB file: {}'.format(args.pfb_out))
	
	print('Making GC model file')
	outdir = os.path.dirname(args.gc_out)
	if not os.path.exists(outdir):
		os.makedirs(outdir)
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
		raise Exception('Generation of GCmodel file with PennCNV failed')
	else:
		print('Created GCmodel file: {}'.format(args.gc_out))

	return 0


def run_snakemake(args):
	
	argv = [
		"-s", os.path.join(SNAKEDIR, "cnv-pipeline.snake"),
		"-p", "-r",
		#"--keep-incomplete",
		"--rerun-incomplete"
	]
	
	argv += [
		'-d', args.directory,
		'--configfile', os.path.join(SNAKEDIR, 'default_config.yaml'), args.config,
		'--config', f'sample_table={args.sample_table}',
			f'snakedir={SNAKEDIR}',
			f'basedir={args.directory}',
			f'configfile={args.config}',
			f'target={args.target}',
	]
	
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
		
	#print(argv)
	
	return snakemake_main(argv)


### Setup ###

def setup_argparse():
	parser = argparse.ArgumentParser(description="run CNV pipeline and helpers")

	group_basic = parser.add_argument_group("General", "General pipeline arguments")

	group_basic.add_argument('--action', '-a', default='run', choices=('run', 'setup-files', 'make-penncnv-files'), help='Action to perform. Default: %(default)s')
	group_basic.add_argument('--config', default='config.yaml', help="Filename of config file. Default: %(default)s")
	group_basic.add_argument('--sample-table', '-s', default='sample_table.txt', help="Filename of sample table. Default: %(default)s")
	#group_basic.add_argument('--data-path', '-p', default='data', help="Filepath to were results are written inside the run directrory. Default: %(default)s")

	group_penncnv = parser.add_argument_group("make-penncnv-files", "Specific arguments for make-penncnv-files")
	group_penncnv.add_argument('--genome', default='GRCh38', choices=('GRCh37', 'GRCh38'),
							   help="Genome build to make the GC model for (uses the files shipped with PennCNV). Default: %(default)s")
	group_penncnv.add_argument('--pfb-out', default='static-data/PennCNV-PFB_from_clusterfile-stats.pfb',
							   help="Filename for generated PFB file. Default: %(default)s")
	group_penncnv.add_argument('--gc-out', default='static-data/PennCNV-GCmodel-GRCh38.gcmodel',
							   help="Filename for generated GCmodel file. Default: %(default)s")

	group_snake = parser.add_argument_group("Snakemake Settings", "Arguments for Snakemake")

	group_snake.add_argument('--target', '-t', default='report', choices=('report', 'cnv-vcf', 'processed-calls', 'PennCNV', 'CBS', 'GADA', 'filtered-data', 'unfiltered-data'),
							 help="Final target of the pipeline. Default: %(default)s")
	group_snake.add_argument('--cluster-profile', '-c', nargs='?', const='cubi-dev', help="Use snakemake profile for job submission to cluster. Default if used: %(const)s")
	group_snake.add_argument('-jobs', '-j', default=20, help="Number of oarallel job submissions in cluster mode. Default: %(default)s")
	group_snake.add_argument('--local-cores', '-n', default=4, help="Number of cores for local submission. Default: %(default)s")
	group_snake.add_argument('--directory', '-d', default=os.getcwd(), help="Directory to run pipeline in. Default: $CWD")
	group_snake.add_argument('snake_options', nargs='*', #argparse.REMAINDER,
							 help="Options to pass to snakemake; separate from normal options with '--'")
	
	return parser
	


if __name__ == '__main__':
	
	args = setup_argparse().parse_args()
	
	#TODO
	# check that conda & PennCNV are set up ?!
	
	if args.action == 'run':
		check_sample_table(args)
		#check_config(args)
		ret = run_snakemake(args)
	elif args.action == 'setup-files':
		ret = copy_setup_files(args)
	elif args.action == 'make-penncnv-files':
		ret = make_penncnv_files(args)

	sys.exit(ret)

