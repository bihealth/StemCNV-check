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
from scripts.py_helpers import *
from scripts.py_exceptions import *

SNAKEDIR = os.path.dirname(os.path.realpath(__file__))
with open(os.path.join(SNAKEDIR, 'default_config.yaml')) as f:
	DEF_CONFIG = yaml.safe_load(f)

### Sanity checks ###
def check_sample_table(args):
	sample_data = read_sample_table(args.sample_table)

	# all_ids = [sid for sid, _, _, _, _ in sample_data]
	# excl = ('no data', 'reference', 'waiting for data')
	# sample_data = [[i, n, p, s, r if r not in excl else ''] for i, n, p, s, r in sample_data] #if all_ids.count(i) == 1

	# Check sample_ids are unique
	#The way the sample_table is read in means that non-unique sample_names will be silently overwrite one another ...
	all_ids = [sid for sid, _, _, _, _ in sample_data]
	non_unique_ids = [sid for sid in all_ids if all_ids.count(sid) > 1]
	if non_unique_ids:
		raise SampleIDNonuniqueError('The following Sample_IDs occur more than once: ' + ', '.join(non_unique_ids))

	# Check sex values
	samples = {sid: sex for sid, _, _, sex, _ in sample_data}
	if any(not s for s in samples.values()):
		missing = [sid for sid, s in samples.items() if not s]
		raise SampleConstraintError("Missing values for 'Sex' in the samplesheet. Affected samples: " + ', '.join(missing))
	elif not all(s in ('m', 'f') for s in map(lambda x: x[0].lower(), samples.values())):
		missing = [f"{sid}: {s}" for sid, s in samples.items() if not s[0].lower() in ('m', 'f')]
		raise SampleConstraintError("Not all values of the 'Sex' column in the samplesheet can be coerced to 'm' or 'f'. Affected samples: " + ', '.join(missing))
	#Check that all reference samples exist
	ref_samples = {rid: sex for _, _, _, sex, rid in sample_data if rid}
	missing_refs = [ref for ref in ref_samples.keys() if ref not in samples.keys()]
	if missing_refs:
		raise SampletableReferenceError("These 'Reference_Sample's do not also exist in the 'Sample_ID' column of the samplesheet: " + ', '.join(missing_refs))
	# Give warning if sex of reference and sample don't match
	sex_mismatch = [f"{s} ({sex})" for s, _, _, sex, ref in sample_data if ref and sex[0].lower() != samples[ref][0].lower()]
	if sex_mismatch:
		sys.stderr.write("Warning: the following samples have a different sex annotation than their Reference_Sample: " + ', '.join(sex_mismatch))
	# Check that Chip_Name & Chip_Pos match the sentrix wildcard regex
	with open(args.config) as f:
		config = yaml.safe_load(f)

	for constraint, val in (('sample_id', 'sid'), ('sentrix_name', 'n'), ('sentrix_pos', 'p')):
		pattern = config_extract(['wildcard_constraints', constraint], config, DEF_CONFIG)
		mismatch = [sid for sid, n, p, _, _ in sample_data if not re.match('^' + pattern + '$', eval(val))]
		if mismatch:
			raise SampleConstraintError(f"The '{constraint}' values for these samples do not fit the expected constraints: " + ', '.join(mismatch))

	sample_data_full = read_sample_table(args.sample_table, with_opt=True)

	# Check optional 'Regions of Interest' column
	if 'Regions_of_Intertest' in sample_data_full[0].keys():
		for data_dict in sample_data_full:
			regions = data_dict['Regions_of_Intertest'].split(';')
			checks = [re.match('^(.*\|)?(chr)?[0-9XY]{1,2}:[0-9]+-[0-9]+$', region) for region in regions]
			if not all(checks):
				raise SampletableRegionError(f"The 'Region_of_Interest' entry for this sample is not properly formatted: {data_dict['Sample_ID']}. Format: (NAME_)?(chr)?[CHR]:[start]-[end], separate multiple regions with only ';'.")


	#Check SNP_clustering extra ids (if defined)
	check_snp_cluster = config_extract(('settings', 'report', 'SNP_comparison', 'include_dendrogram'), config, DEF_CONFIG)
	extra_sample_def = config_extract(('settings', 'report', 'SNP_comparison', 'extra_samples'), config, DEF_CONFIG)
	if check_snp_cluster:
		for sample_id in samples.keys():
			cluster_ids = collect_SNP_cluster_ids(sample_id, extra_sample_def, sample_data_full)
			missing = [sid for sid in cluster_ids if sid not in samples.keys()]
			if missing:
				raise SampletableReferenceError(
					"These Sample_ID defined as controls for SNP clustering dendrogram (or derived from one of the columns there) do not exist in the 'Sample_ID' column of the samplesheet: " + ', '.join(
						missing))


# Assume that the default_config hasn't been altered & is correct
def check_config(args):

	with open(args.config) as f:
		config = yaml.safe_load(f)

	if config_extract(['settings', 'make_cnv_vcf', 'mode'], config, DEF_CONFIG) not in ('combined-calls', 'split-tools'):
		raise ConfigValueError(
			'Value not allowed for settings$make.cnv.vcf$mode: "{}"'.format(config['settings']['make.cnv.vcf']['mode']))

	#TODO:
	# - check that required values are filled in
	# - check that value types are correct / match config (incl list len?)
	# - maybe: us the same apprach that seasnap has for this?


def check_installation():
	pass


def make_PennCNV_sexfile(args):
	"""Make the `sexfile` that PennCNV requirs from the sampletable & config. Needs to updated for each run,
	but would trigger snakemake reruns if done as a snakemake rule"""
	#Description of needed format:
	# A 2-column file containing filename and sex (male/female) for sex chromosome calling with -chrx argument. The first
	# tab-delimited column should be the input signal file name, while the second tab-delimited column should be male or female.
	# Alternatively, abbreviations including m (male), f (female), 1 (male) or 2 (female) are also fine.

	with open(args.config) as f:
		config = yaml.safe_load(f)

	basepath = args.directory
	datapath = config['data_path']
	outfilename = os.path.join(basepath, "penncnv-sexfile.txt")

	filter = config_extract(['settings', 'PennCNV', 'filter-settings'], config, DEF_CONFIG)

	sample_data = read_sample_table(args.sample_table)

	with open(outfilename, 'w') as f:
		for sample_id, _, _, sex, _ in sample_data:
			inputfile = os.path.join(basepath, datapath, f"{sample_id}", f"{sample_id}.filtered-data.{filter}.tsv")
			#ensure its consistently 'm'/'f'
			sex = sex.lower()[0]
			f.write(f"{inputfile}\t{sex}\n")



### Actions ###

def copy_setup_files(args):
	#TODO: copy into PWD or into directroy?
	#TODO: add different levels of copying config (only required, `normal` and full)
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
	
	vcf_files = [os.path.join(args.directory, datapath, f"{sample_id}", f"{sample_id}.unprocessed.vcf") for sample_id, _, _, _, _ in sample_data]
	vcf_present = [vcf for vcf in vcf_files if os.path.exists(vcf)]
	
	if vcf_present:
		use_vcf = vcf_present[0]
	else:
		use_vcf = vcf_files[0]
		print('Running snakemake to get:', use_vcf)
		args.snake_options += [	use_vcf ]
		
		#TODO: The code terminate here because snakemake.main includes a sys.exit() call
		# --> need to use snakemake.snakemake and parse args (manually?) instead
		ret = run_snakemake(args)
		
		
# 		snakemake(os.path.join(SNAKEDIR, "cnv-pipeline.snake"), 
# 		
# 		            cores=args.cores,
#             local_cores=args.local_cores,
#             nodes=args.jobs,
#             
#             config=config,
#             configfiles=args.configfile,
#             config_args=args.config,
#             workdir=args.directory,
#             
#             		'--configfile', os.path.join(SNAKEDIR, 'default_config.yaml'), args.config,
# 		'--config', f'sample_table={args.sample_table}',
# 			f'snakedir={SNAKEDIR}',
# 			f'basedir={args.directory}',
# 			f'configfile={args.config}',
# 			f'target={args.target}',
#             
#             
#             targets= use_vcf,
#             )
		
		if ret:
			raise Exception('Snakemake run to get vcf failed')
	
	#TODO logger calls
	print(f'Making PFB file: {args.pfb_out}')
	outdir = os.path.dirname(args.pfb_out)
	if not os.path.exists(outdir):
		os.makedirs(outdir)

	argv = ['Rscript', f"{SNAKEDIR}/scripts/make_PFB_from_vcf.R", use_vcf, args.pfb_out]
	print( 'Executing: `'  + ' '.join(argv) + '`')
	ret = subprocess.call(argv)
	
	if ret:
		raise Exception('Generation of PFB from vcf file failed')
	else:
		print('Created PFB file: {}'.format(args.pfb_out))
	
	print(f'Making GC model file: {args.gc_out}')
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
	print( 'Executing: `'  + ' '.join(argv) + '`')
	ret = subprocess.call(argv)

	if ret:
		raise Exception('Generation of GCmodel file with PennCNV failed')
	else:
		print('Created GCmodel file: {}'.format(args.gc_out))

	return 0


def run_snakemake(args):

	make_PennCNV_sexfile(args)
	
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
	group_basic.add_argument('--config', '-c', default='config.yaml', help="Filename of config file. Default: %(default)s")
	group_basic.add_argument('--sample-table', '-s', default='sample_table.txt', help="Filename of sample table. Default: %(default)s")
	#group_basic.add_argument('--data-path', '-p', default='data', help="Filepath to were results are written inside the run directrory. Default: %(default)s")

	group_penncnv = parser.add_argument_group("make-penncnv-files", "Specific arguments for make-penncnv-files")
	group_penncnv.add_argument('--genome', default='GRCh38', choices=('GRCh37', 'GRCh38'),
							   help="Genome build to make the GC model for (uses the files shipped with PennCNV). Default: %(default)s")
	group_penncnv.add_argument('--pfb-out', default='static-data/PennCNV-PFB_from_clusterfile-stats_{genome}.pfb',
							   help="Filename for generated PFB file. Default: %(default)s")
	group_penncnv.add_argument('--gc-out', default='static-data/PennCNV-GCmodel-{genome}.gcmodel',
							   help="Filename for generated GCmodel file. Default: %(default)s")

	group_snake = parser.add_argument_group("Snakemake Settings", "Arguments for Snakemake")

	group_snake.add_argument('--target', '-t', default='report', choices=('report', 'cnv-vcf', 'processed-calls', 'PennCNV', 'CBS', 'GADA', 'SNP-probe-data'),
							 help="Final target of the pipeline. Default: %(default)s")
	group_snake.add_argument('--cluster-profile', '-p', nargs='?', const='cubi-dev', help="Use snakemake profile for job submission to cluster. Default if used: %(const)s")
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
		args.pfb_out = args.pfb_out.format(genome=args.genome)
		args.gc_out = args.gc_out.format(genome=args.genome)
		ret = make_penncnv_files(args)

	sys.exit(ret)

