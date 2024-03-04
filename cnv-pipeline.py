#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Wrapper around snakemake to kick-off pipeline
"""

import argparse
from loguru import logger

import os
import re
import subprocess
import sys
import tempfile
import yaml
from collections import defaultdict
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
	# Check that all reference samples exist
	ref_samples = {rid: sex for _, _, _, sex, rid in sample_data if rid}
	missing_refs = [ref for ref in ref_samples.keys() if ref not in samples.keys()]
	if missing_refs:
		raise SampletableReferenceError("These 'Reference_Sample's do not also exist in the 'Sample_ID' column of the samplesheet: " + ', '.join(missing_refs))
	# Give warning if sex of reference and sample don't match
	sex_mismatch = [f"{s} ({sex})" for s, _, _, sex, ref in sample_data if ref and sex[0].lower() != samples[ref][0].lower()]
	if sex_mismatch:
		raise SampletableReferenceError("The following samples have a different sex annotation than their Reference_Sample: " + ', '.join(sex_mismatch))
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


#Note: this is WIP
#interactive prompts to check missing files (pfb, gcmodel, genome info, etc), then create them
def interactive_check_staticdata(args, config):

	for req in ('bpm_manifest_file', 'egt_cluster_file'):
		if not req in config['static_data'] or not config['static_data'][req]:
			raise InputFileError(f"Required config entry is missing: static-data:{req}")
		if not os.path.isfile(config['static_data'][req]):
			raise InputFileError(f"Defined static data file for '{req}' does not exist.")

	for req in ('csv_manifest_file', 'genome_fasta_file', 'genome_gtf_file'):
		if not req in config['static_data'] or not config['static_data'][req]:
			logger.warning(f"No static input for '{req}' is defined. Generated vcf file will be missing some features.")
		if req in config['static_data'] and not os.path.isfile(config['static_data'][req]):
			raise InputFileError(f"Defined static data file for '{req}' does not exist.")

	missing_creatable = set()
	for req in ('csv_manifest_file', 'genome_fasta_file', 'genome_gtf_file'):
		if not req in config['static_data'] or not config['static_data'][req]:
			missing_creatable.add(req)

		if req in config['static_data'] and not os.path.isfile(config['static_data'][req]):
			missing_creatable.add(req)
			logger.warning(f"Defined static data file for '{req}' does not exist. This file can be created")

	if missing_creatable:
		logger.warning("The following required static data files are not defined in the config or do not exist yet, but they can be auto-generated: {}".format(', '.join(missing_creatable)))
		ask = input("Should the missing, auto-generatable files be created now? [yN]")
		if not ask or ask[0].lower() != 'y':
			logger.warning("Not creating missing required files, can not continue.")
			raise InputFileError("Required static data files are missing: " + ', '.join(missing_creatable))

		genome = args.genome

	rewrite_config = False


def check_config(args):

	with open(args.config) as f:
		config = yaml.safe_load(f)
	with open(os.path.join(SNAKEDIR, 'allowedvalues_config.yaml')) as f:
		allowed_values = yaml.safe_load(f)

	## Check required enrties

	# Files: static-data/*
	for req in allowed_values['static_data'].keys():
		if req in ('pfb_file', 'GCmodel_file', 'genomeInfo_file', 'array_density_file', 'array_gaps_file'):
			infostr = " You can create it by running `cnv-pipeline -a make-staticdata`"
		else:
			infostr = ""
		if not req in config['static_data'] or not config['static_data'][req]:
			raise InputFileError(f"Required config entry is missing: static-data:{req}")
		if not os.path.isfile(config['static_data'][req]):
			raise InputFileError(f"Static data file '{req}' does not exist." + infostr)

	# Folders: log, data, raw-input
	for req in ('data_path', 'log_path', 'raw_data_folder'):
		try:
			if not config[req]:
				raise ConfigValueError(f"Required config entry is missing: {req}")
			if not os.path.isdir(config[req]):
				warn_str = f"Required Entry '{req}' is not an existing folder! Attempting to create it."
				logger.warning(warn_str, ConfigValueWarning)
				os.makedirs(config[req], exist_ok=False)
		except KeyError:
			raise ConfigValueError(f"Required config entry is missing: {req}")
	# Other settings: reports/*/filetype
	if not 'reports' in config:
		logger.warning('No reports are defined in the config, only tabular & vcf files will be created!', ConfigValueWarning)
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

	defined_filtersets = set(config_extract(['settings', 'probe-filter-sets'], config, DEF_CONFIG).keys())
	#defined_filtersets = set(config['settings']['probe-filter-sets'].keys()) | set(DEF_CONFIG['settings']['probe-filter-sets'].keys())
	check_functions = defaultdict(lambda: lambda x, v: bool(re.match('^'+str(v)+'$', str(x))))
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
			logger.warning(warn_str, ConfigValueWarning)
		raise ConfigValueError('The config contains values that are not allowed')


## Helper Functions ##

def make_singularity_args(config, not_existing_ok=False):
	"""Collect all outside filepaths that need to be bound inside container"""

	bind_points = [
		(config['data_path'], '/outside/data'),
		(config['raw_data_folder'], '/outside/rawdata'),
		(config['log_path'], '/outside/logs'),
		(SNAKEDIR, '/outside/snakedir')
	]

	for name, file in config['static_data'].items():
		# Can only mount existing files
		if not_existing_ok and not os.path.isfile(file):
			continue
		bind_points.append((file, '/outside/static/{}'.format(os.path.basename(file))))

	return "-B " + ','.join(f"{host}:{cont}" for host, cont in bind_points)


# This is done outside of snakemake so that the file can be updated without this triggering reruns
# (which it would independently of mtime if created by a rule)
def make_PennCNV_sexfile(args):
	"""Make the `sexfile` that PennCNV requires from the sampletable & config. Needs to be updated for each run,
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
	if filter == '__default__':
		filter = config_extract(['settings', 'default-filter-set'], config, DEF_CONFIG)

	sample_data = read_sample_table(args.sample_table)

	with open(outfilename, 'w') as f:
		for sample_id, _, _, sex, _ in sample_data:
			inputfile = os.path.join(basepath, datapath, f"{sample_id}", f"{sample_id}.filtered-data.{filter}.tsv")
			# Ensure its consistently 'm'/'f'
			sex = sex.lower()[0]
			f.write(f"{inputfile}\t{sex}\n")

### Actions ###

def copy_setup_files(args):
	subprocess.call(['cp', f"{SNAKEDIR}/default_config.yaml", args.config])
	subprocess.call(['cp', f"{SNAKEDIR}/sample_table.txt", args.sample_table])
	logger.info('Created example files: ...')


def create_missing_staticdata(args):
	# Check if any vcf file is present, generate one if none are
	# This is done because the vcf files generated with gtc2vcf contain the full information
	# about GT (background) counts from the egt clusterfile, that can be used to make the pfb file for PennCNV
	sample_data = read_sample_table(args.sample_table)
	with open(args.config) as f:
		config = yaml.safe_load(f)
	datapath = config_extract(('data_path',), config, DEF_CONFIG)

	vcf_files = [os.path.join(args.directory, datapath, f"{sample_id}", f"{sample_id}.unprocessed.vcf") for
				 sample_id, _, _, _, _ in sample_data]
	vcf_present = [vcf for vcf in vcf_files if os.path.exists(vcf)]

	if vcf_present:
		use_vcf = vcf_present[0]
	else:
		use_vcf = vcf_files[0]
		print('Running snakemake to get:', use_vcf)

		ret = snakemake(
			os.path.join(SNAKEDIR, "cnv-pipeline.smk"),
			local_cores=args.local_cores,
			cores=args.local_cores,
			workdir=args.directory,
			configfiles=[os.path.join(SNAKEDIR, 'default_config.yaml'), args.config],
			printshellcmds=True,
			force_incomplete=True,
			use_conda=True,
			use_singularity=not args.no_singularity,
			singularity_args='' if args.no_singularity else make_singularity_args(config, True),
			config={'sample_table': args.sample_table,
					'snakedir': SNAKEDIR,
					'basedir': args.directory,
					'configfile': args.config,
					'use_singularity': not args.no_singularity
					},
			targets=[use_vcf]
		)

		if not ret:
			logger.error('Snakemake run to get vcf failed')
			sys.exit(1)

	# Run extra snakemake to check which files are missing & create them accordingly
	with tempfile.TemporaryDirectory() as tmpdir:
		density_windows = config_extract(('settings', 'postprocessing', 'density.windows',), config, DEF_CONFIG)
		min_gap_size = config_extract(('settings', 'postprocessing', 'min.gap.size',), config, DEF_CONFIG)
		ret = snakemake(
			os.path.join(SNAKEDIR, "staticdata_creation.smk"),
			local_cores=args.local_cores,
			cores=args.local_cores,
			workdir=args.directory,
			use_singularity=not args.no_singularity,
			use_conda=True,
			conda_frontend=args.conda_frontend,
			printshellcmds=True,
			force_incomplete=True,
			config={'snakedir': SNAKEDIR,
					'use_singularity': not args.no_singularity,
					'TMPDIR': tmpdir,
					'genome': args.genome,
					'vcf_input_file': use_vcf,
					'pfb_outname': args.penncnv_pfb_out,
					'gcmodel_outname': args.penncnv_gc_out,
					'chrominfo_outname': args.chromosome_info_out,
					'array_density_outname': args.array_density_out,
					'array_gaps_outname': args.array_gaps_out,
					'density_windows': density_windows,
					'min_gap_size': min_gap_size,
					}
		)

	logger.info("""All files generated, add the following lines to the static-data section of your config file:
  pfb_file: {pfb_out}
  GCmodel_file: {gcmodel_out}
  array_density: {density_out}
  array_gaps: {gap_out}
  genomeInfo_file: {info_out}""".format(
		pfb_out=args.penncnv_pfb_out, gcmodel_out=args.penncnv_gc_out, density_out=args.array_density_out,
		gap_out=args.array_gaps_out, info_out=args.chromosome_info_out)
	)

	return ret


def run_snakemake(args):

	# Ensure that sexfile for PennCNV exists
	make_PennCNV_sexfile(args)
	with open(args.config) as f:
		config = yaml.safe_load(f)
	
	argv = [
		"-s", os.path.join(SNAKEDIR, "cnv-pipeline.smk"),
		"-p", "--rerun-incomplete",
		"--use-conda", "--conda-frontend", args.conda_frontend,
	]
	if args.no_singularity:
		logger.warning("Running without singularity containers is a legacy feature, pipeline will fail without local PennCNV installation based on the install.sh script!")
		use_singularity = False
	else:
		argv += [ "--use-singularity",
				  "--singularity-args", make_singularity_args(config)
				  ]
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
	
	return main(argv)


### Setup ###

def setup_argparse():
	parser = argparse.ArgumentParser(description="run CNV pipeline and helpers")

	group_basic = parser.add_argument_group("General", "General pipeline arguments")

	group_basic.add_argument('--action', '-a', default='run', choices=('run', 'setup-files', 'make-staticdata'), help='Action to perform. Default: %(default)s')
	group_basic.add_argument('--config', '-c', default='config.yaml', help="Filename of config file. Default: %(default)s")
	group_basic.add_argument('--sample-table', '-s', default='sample_table.txt', help="Filename of sample table. Default: %(default)s")
	group_basic.add_argument('--directory', '-d', default=os.getcwd(),
							 help="Directory to run pipeline in. Default: $CWD")

	# group_basic.add_argument('--non-interactive', '-y', action="store_true",
	# 						 help="Skip all interactive questions of the wrapper")
	group_basic.add_argument('--conda-frontend', default='mamba', choices=('mamba', 'conda'), help="Conda frontend to use. Default: %(default)s")
	group_basic.add_argument('--no-singularity', action='store_true',
							 help="Do not use singularity/docker, you will need a local PennCNV installation instead (see install.sh)")
	group_basic.add_argument('--verbose', '-v', action='store_true', help="Verbose output")

	group_static = parser.add_argument_group("make-staticdata", "Details and file naming for make-staticdata")
	group_static.add_argument('--genome', default='hg38', choices=('hg19', 'hg38'), help="Genome build to use (UCSC names). Default: %(default)s")
	group_static.add_argument('--snp-array-name', default=None, help="A name or identifier string for the snp-array, can used in filesnames. No Default.")
	group_static.add_argument('--penncnv-pfb-out', default='static-data/PennCNV-PFB_{genome}{array}.pfb',
							   help="Filename for generated PFB file. Default: %(default)s")
	group_static.add_argument('--penncnv-gc-out', default='static-data/PennCNV-GCmodel_{genome}{array}.gcmodel',
							   help="Filename for generated GCmodel file. Default: %(default)s")
	group_static.add_argument('--array-density-out', default='static-data/density_{genome}{array}.bed',
							   help="Filename for generated bed file with probe density. Default: %(default)s")
	group_static.add_argument('--array-gaps-out', default='static-data/gaps_{genome}{array}.bed',
							  help="Filename for generated bed file with probe gaps. Default: %(default)s")
	group_static.add_argument('--chromosome-info-out', default='static-data/UCSC_{genome}_chromosome-info.tsv',
							   help="Filename for generated chromosome info file. Default: %(default)s")

	group_snake = parser.add_argument_group("Snakemake Settings", "Arguments for Snakemake (also affects make-staticdata)")

	group_snake.add_argument('--target', '-t', default='report', choices=('report', 'cnv-vcf', 'processed-calls', 'PennCNV', 'CBS', 'SNP-probe-data'),
							 help="Final target of the pipeline. Default: %(default)s")
	group_snake.add_argument('--cluster-profile', '-p', nargs='?', const='cubi-dev', help="Use snakemake profile for job submission to cluster. Default if used: %(const)s")
	group_snake.add_argument('-jobs', '-j', default=20, help="Number of oarallel job submissions in cluster mode. Default: %(default)s")
	group_snake.add_argument('--local-cores', '-n', default=4, help="Number of cores for local submission. Default: %(default)s")
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
	elif args.action == 'make-staticdata':
		args.snp_array_name = '_' + args.snp_array_name if args.snp_array_name and not args.snp_array_name[0] in ".-_" else ""
		args.penncnv_pfb_out = args.penncnv_pfb_out.format(genome=args.genome, array=args.snp_array_name)
		args.penncnv_gc_out = args.penncnv_gc_out.format(genome=args.genome, array=args.snp_array_name)
		args.array_density_out = args.array_density_out.format(genome=args.genome, array=args.snp_array_name)
		args.array_gaps_out = args.array_gaps_out.format(genome=args.genome, array=args.snp_array_name)
		args.chromosome_info_out = args.chromosome_info_out.format(genome=args.genome)
		if not os.path.isdir(args.directory):
			os.makedirs(args.directory)
		ret = create_missing_staticdata(args)

	sys.exit(ret)

