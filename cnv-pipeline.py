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
import warnings
import yaml
from collections import defaultdict
# import ruamel.yaml as ruamel_yaml
# from snakemake import RERUN_TRIGGERS
from snakemake import main, snakemake
from scripts.py_helpers import *
from scripts.py_exceptions import *

SNAKEDIR = os.path.dirname(os.path.realpath(__file__))
with open(os.path.join(SNAKEDIR, 'default_config.yaml')) as f:
	DEF_CONFIG = yaml.safe_load(f)

### Sanity checks ###
def check_sample_table(args):
	sample_data = read_sample_table(args.sample_table)

	# Check sample_ids are unique
	#The way the sample_table is read in means that non-unique sample_ids would silently overwrite one another ...
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
		warnings.warn("The following samples have a different sex annotation than their Reference_Sample: " + ', '.join(sex_mismatch))
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


# Assume that the default_config hasn't been altered & is correct
def check_config(args):

	with open(args.config) as f:
		config = yaml.safe_load(f)
	with open(os.path.join(SNAKEDIR, 'allowedvalues_config.yaml')) as f:
		allowed_values = yaml.safe_load(f)

	## Check required enrties
	# Folders: log, data, raw-input
	for req in ('data_path', 'log_path', 'raw_data_folder'):
		try:
			if not config[req]:
				raise ConfigValueError(f"Required config entry is missing: {req}")
			if not os.path.isdir(config[req]):
				warn_str = f"Required Entry '{req}' is not an existing folder! Attempting to create it."
				warnings.warn(warn_str, ConfigValueWarning)
				os.makedirs(config[req], exist_ok=False)
		except KeyError:
			raise ConfigValueError(f"Required config entry is missing: {req}")
	# Files: static-data/*
	for req in allowed_values['static_data'].keys():
		try:
			if not config['static_data'][req]:
				raise ConfigValueError(f"Required config entry is missing: static-data:{req}")
			if not os.path.isfile(config['static_data'][req]):
				if req in ('pfb_file', 'GCmodel_file'):
					info = " You can create it by running `cnv-pipeline -a make-penncnv-files`"
				else:
					info = ""
				raise ConfigValueError(f"Required static data file '{req}' is missing." + info)
		except KeyError:
			raise ConfigValueError(f"Required config entry is missing: static-data:{req}")
	# Other settings: reports/*/filetype
	if not 'reports' in config:
		#TODO: test if the pipeline won't crash in this case
		warnings.warn('No reports are defined in the config, only tabular & vcf files will be created!', ConfigValueWarning)
	else:
		for rep in config['reports']:
			if rep == '__default__':
				continue
			try:
				if not config['reports'][rep]['file_type']:
					raise ConfigValueError(f"Required config entry is missing: reports:{rep}:file_type")
			except KeyError:
				raise ConfigValueError(f"Required config entry is missing: reports:{rep}:file_type")

	## Check value ranges
	sample_data = read_sample_table(args.sample_table, with_opt=True)

	def parse_scientific(n):
		if isinstance(n, (int, float)):
			return n
		if re.match('^[0-9.]+e[0-9]+$', n):
			n, e = n.split('e')
			out = float(n) * 10**int(e)
			return int(out) if int(out) == out else out
		else:
			ValueError(f"{n} can not be coerced to a number")

	defined_filtersets = set(config['settings']['probe-filter-sets'].keys()) | set(DEF_CONFIG['settings']['probe-filter-sets'].keys())
	check_functions = defaultdict(lambda: lambda x, v: bool(re.match('^'+str(v)+'$', x)))
	def_functions = {
		'list': lambda x, v: type(x) == list,
		'str': lambda x, v: type(x) == str,
		'int': lambda x, v: type(parse_scientific(x)) == int,
		'float': lambda x, v: type(parse_scientific(x)) == float or type(parse_scientific(x)) == int,
		'bool': lambda x, v: type(x) == bool,
		'len': lambda x, v: len(x) == v,
		'le': lambda x, v: parse_scientific(x) <= v,
		'ge': lambda x, v: parse_scientific(x) >= v,
		'filterset': lambda x, v: x in defined_filtersets or x == '__default__',
		'filtersetnodefault': lambda x, v: x in defined_filtersets,
		'sections': lambda x, v: x in allowed_values['allowed_sections'],
		'sectionsall': lambda x, v: x == '__all__' or all(i in allowed_values['allowed_sections'] for i in x),
		'insamplesheet': lambda x, v: re.sub('^_+', '', x) in sample_data[0].keys()
	}
	check_functions.update(def_functions)

	#Helper function:
	def formatfunc(s):
		strpart = s.rstrip('0123456789')
		if strpart in def_functions:
			numpart = int(s[len(strpart):]) if s[len(strpart):] else None
		else: # regex function
			numpart = strpart
		return strpart, numpart

	# Check all config entries
	errors = []
	for flatkey, config_value in flatten(config).items():
		# Need to change the key for variable config sections
		flatkey_ = re.sub('reports:[^:]+', 'reports:__report', flatkey)
		flatkey_ = re.sub('tools:[^:]+', 'tools:__tool', flatkey_)
		flatkey_ = re.sub('settings:probe-filter-sets:[^:]+', 'settings:probe-filter-sets:__filterset', flatkey_)
		funcs = config_extract(flatkey_.split(':'), allowed_values, allowed_values)
		if funcs is None:
			print(flatkey, config_value)
			continue
		funcs, *list_funcs = funcs.split('__')
		for func_key in funcs.split('_'):
			func_key, func_value = formatfunc(func_key)
			check = check_functions[func_key](config_value, func_value)
			if not check:
				errors.append((flatkey, config_value, func_key, func_value))
		for func_key in list_funcs:
			func_key, func_value = formatfunc(func_key)
			check = all(check_functions[func_key](i, func_value) for i in config_value)
			if not check:
				errors.append((flatkey, config_value, func_key, func_value))

	#TODO: better help message where regex is a list of values
	#help_strings = defaultdict(lambda: lambda v: "only these entries: " + v[1:-1].replace('|', ', '))
	help_strings = defaultdict(lambda: lambda v: "matching this regex: " + v)
	help_strings.update({
		'list': lambda v: 'a list with ',
		'str': lambda v: 'characters (add quotes for numbers or similar)',
		'int': lambda v: 'integers (whole numbers)',
		'float': lambda v: 'numbers',
		'bool': lambda v: 'booleans (True/False)',
		'len': lambda v: f"exactly {v} entries",
		'le': lambda v: f"numbers <={v}",
		'ge': lambda v: f"numbers >={v}",
		'filterset': lambda v: "the defined filtersets ({}) or __default__".format(', '.join(defined_filtersets)),
		'filtersetnodefault': lambda v: "only the defined filtersets ({}) but not '__default__'".format(', '.join(defined_filtersets)),
		'sections': lambda v: "the defined report sections: " + ', '.join(allowed_values['allowed_sections']),
		'sectionsall': lambda v: "either '__all__' or a list with of defined report sections: " + ', '.join(allowed_values['allowed_sections']),
		'insamplesheet': lambda v: "column names of the samplesheet: " + ', '.join(sample_data[0].keys())
	})
	if errors:
		for flatkey, config_value, func_key, func_value in errors:
			warn_str = f"The config entry '{config_value}' for '{flatkey}' is invalid. Allowed value" + \
					   ('s are ' if func_key != 'list' else ' is a list with ') + \
					   help_strings[func_key](func_value) + '.'
			warnings.warn(warn_str, ConfigValueWarning)
		raise ConfigValueError('The config contains values that are not allowed')


def check_installation():
	pass


def make_PennCNV_sexfile(args):
	"""Make the `sexfile` that PennCNV requires from the sampletable & config. Needs to be updated for each run,
	but would trigger snakemake reruns if done as a snakemake rule"""
	#TODO: would a subworkflow still trigger reruns?
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
	if filter == '__default__':
		filter = config_extract(['settings', 'default-filter-set'], config, DEF_CONFIG)

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

#TODO: maybe it'd be better if this was part of the actual snakemake workflow?
# -> the issue though is that it shouldn't depend on any given sample from the view of snakemake
def make_penncnv_files(args):

	# Check if any vcf file is present
	sample_data = read_sample_table(args.sample_table)
	with open(args.config) as f:
		config = yaml.safe_load(f)
	datapath = config_extract(('data_path', ), config, DEF_CONFIG)
	
	vcf_files = [os.path.join(args.directory, datapath, f"{sample_id}", f"{sample_id}.unprocessed.vcf") for sample_id, _, _, _, _ in sample_data]
	vcf_present = [vcf for vcf in vcf_files if os.path.exists(vcf)]
	
	if vcf_present:
		use_vcf = vcf_present[0]
	else:
		use_vcf = vcf_files[0]
		print('Running snakemake to get:', use_vcf)
		#args.snake_options += [	use_vcf ]

		ret = snakemake(
			os.path.join(SNAKEDIR, "cnv-pipeline.smk"),
			local_cores=args.local_cores,
			cores=args.local_cores,
			workdir=args.directory,
			configfiles=[os.path.join(SNAKEDIR, 'default_config.yaml'), args.config],
			printshellcmds=True,
			force_incomplete=True,
			config={'sample_table': args.sample_table,
					'snakedir': SNAKEDIR,
					'basedir': args.directory,
					'configfile': args.config
					},
			targets=[use_vcf]
		)

		if not ret:
			raise Exception('Snakemake run to get vcf failed')


	#TODO: externalise this to an additional snakefile? (-> easiest way to use the pennCNV docker container?)

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

	# Ensure that sexfile for PennCNV exists
	make_PennCNV_sexfile(args)
	
	argv = [
		"-s", os.path.join(SNAKEDIR, "cnv-pipeline.smk"),
		"-p", #"-r", is default now
		"--rerun-incomplete"
	]
	if args.no_singularity:
		warnings.warn("Usage of singularity/docker by snakemake is disabled, pipeline will fail without local PennCNV installation based on the install.sh script!")
		use_singularity = False
	else:
		argv += [ "--use-singularity" ]
		use_singularity = True
	
	argv += [
		'-d', args.directory,
		'--configfile', os.path.join(SNAKEDIR, 'default_config.yaml'), args.config,
		'--config', f'sample_table={args.sample_table}',
			f'snakedir={SNAKEDIR}',
			f'basedir={args.directory}',
			f'configfile={args.config}',
			f'target={args.target}',
			f'use_singularity={use_singularity}'
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
	
	return main(argv)


### Setup ###

def setup_argparse():
	parser = argparse.ArgumentParser(description="run CNV pipeline and helpers")

	group_basic = parser.add_argument_group("General", "General pipeline arguments")

	group_basic.add_argument('--action', '-a', default='run', choices=('run', 'setup-files', 'make-penncnv-files'), help='Action to perform. Default: %(default)s')
	group_basic.add_argument('--config', '-c', default='config.yaml', help="Filename of config file. Default: %(default)s")
	group_basic.add_argument('--sample-table', '-s', default='sample_table.txt', help="Filename of sample table. Default: %(default)s")

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
	group_snake.add_argument('--no-singularity', action='store_true', help="Do not use singularity/docker, you will need a local PennCNV installation instead (see install.sh)")
	group_snake.add_argument('snake_options', nargs='*', #argparse.REMAINDER,
							 help="Options to pass to snakemake; separate from normal options with '--'")
	
	return parser
	


if __name__ == '__main__':
	
	args = setup_argparse().parse_args()

	if args.action == 'run':
		check_sample_table(args)
		check_config(args)
		if not os.path.isdir(args.directory):
			os.makedirs(args.directory)
		ret = run_snakemake(args)
	elif args.action == 'setup-files':
		ret = copy_setup_files(args)
	elif args.action == 'make-penncnv-files':
		args.pfb_out = args.pfb_out.format(genome=args.genome)
		args.gc_out = args.gc_out.format(genome=args.genome)
		if not os.path.isdir(args.directory):
			os.makedirs(args.directory)
		ret = make_penncnv_files(args)

	sys.exit(ret)

