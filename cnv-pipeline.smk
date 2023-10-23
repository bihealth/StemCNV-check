# -*- coding: utf-8 -*-

import os
import tempfile
import yaml
from scripts.py_helpers import *
from scripts.py_exceptions import *
# Configuration ================================================================


if not config:
  # This should not happen unless the snakefile is invoked directly from cmd-line
  # try to load a baseline config file
  CONFIGFILE = "config.yaml"
  configfile: CONFIGFILE
  removetempconfig = False
else:
  # Save config to tempfile & use that instead, this ensures the combined default + user configs are used
  # It also archives all options passed by command line and allows them all to be displayed in report
  f, CONFIGFILE = tempfile.mkstemp(suffix = '.yaml', text=True)
  os.close(f)
  removetempconfig = True
  with open(CONFIGFILE, 'w') as yamlout:
    yaml.dump(config, yamlout)


SAMPLETABLE = config['sample_table'] if 'sample_table' in config else 'sample_table.txt' # Defined by wrapper
SNAKEDIR = config['snakedir'] if 'snakedir' in config else os.path.dirname(os.path.realpath(__file__)) #Defined by wrapper
BASEPATH = config['basedir'] if 'basedir' in config else os.getcwd() #Defined by wrapper
DATAPATH = config['data_path'] if os.path.isabs(config['data_path']) else os.path.join(BASEPATH, config['data_path'])
LOGPATH = config['log_path'] if os.path.isabs(config['log_path']) else os.path.join(BASEPATH, config['log_path'])
TARGET = config['target'] if 'target' in config else 'report' #Defined by wrapper
IDAT_INPUT = config['raw_data_folder']

wildcard_constraints:
  # This seems to cause issues -> will also make debugging harder for users
  #filter=list(config['settings']['probe-filter-sets'].keys()),
  sample_id=config['wildcard_constraints']['sample_id'],
  #sentrix_name=config['wildcard_constraints']['sentrix_name'],
  #sentrix_pos=config['wildcard_constraints']['sentrix_pos'],

#Never submit these to cluster
localrules:
  relink_gencall,
  all

sample_data = read_sample_table(SAMPLETABLE)
sample_data_full = read_sample_table(SAMPLETABLE, with_opt=True)
#print(sample_data)

# Helper functions =============================================================

def get_tool_filter_settings(tool):
  if tool.split(':')[0] == 'report':
    report_settings = config['reports'][tool.split(':')[1]]
    out = config_extract((tool.split(':')[2], 'filter-settings'), report_settings, config['reports']['__default__'])
  else:
    out = config['settings'][tool]['filter-settings']
  if out == '__default__':
    out = config['settings']['default-filter-set']
  return out


def get_tool_resource(tool ,resource):
  if not resource in ('threads', 'memory', 'runtime', 'partition', 'cmd-line-params'):
    raise KeyError(f"This resource can not be defined: {resource}")
  if not tool in config['tools']:
    return config['tools']['__default__'][resource]
  else:
    if resource in config['tools'][tool]:
      return config['tools'][tool][resource]
    else:
      return config['tools']['__default__'][resource]

# Rules ========================================================================

def get_cnv_vcf_output(mode):
  addition = config['settings']['make_cnv_vcf']['name_addition']
  if mode == 'combined-calls':
    return os.path.join(DATAPATH, "{sample_id}", f"{{sample_id}}.combined-cnv-calls{addition}.vcf")
  elif mode == 'split-tools':
    return [os.path.join(DATAPATH, "{sample_id}", f"{{sample_id}}.{tool}-cnv-calls{addition}.vcf")
            for tool in config['settings']['CNV.calling.tools']]
  else:
    raise ConfigValueError('Value not allowed for settings$make.cnv.vcf$mode: "{}"'.format(config['settings']['make_cnv_vcf']['mode']))

def get_target_files():
  #TODO/maybe:
  # allow multiple targets to be specified?
  #
  #Target options: ('report', 'cnv-vcf', 'processed-calls', 'PennCNV', 'CBS', 'GADA', 'SNP-probe-data'),
  all_samples = [sample_id for sample_id, _, _, _, _ in sample_data]
  # Report
  if TARGET == 'report':
    return expand(os.path.join(DATAPATH, "{sample_id}", "{sample_id}.{report_filetype}"),
                  sample_id = all_samples,
                  report_filetype = [rep+'.'+config['reports'][rep]['file_type'] for rep in config['reports'].keys() if rep != '__default__']) #+ \
           # expand(os.path.join(DATAPATH,"{sample_id}","{sample_id}.summary-check.tsv"),
           #        sample_id = all_samples)
  # cnv-vcf
  if TARGET == 'cnv-vcf':
    return expand(get_cnv_vcf_output(config['settings']['make_cnv_vcf']['mode']),
                  sample_id = all_samples)
      #      )
      #
      # [os.path.join(BASEPATH, DATAPATH, f"{sample_id}", f"{sample_id}.{{tool}}-cnv-calls.{filter}.vcf") for
      #              sample_id, _, _, _, _ in sample_data],
      #             tool = 'combined' if config['settings']['make_cnv_vcf']['mode'] == 'combined-calls' else config['settings']['CNV.calling.tools'])

  # Target Processed-calls
  if TARGET == 'processed-calls':
    return expand(os.path.join(DATAPATH, "{sample_id}", "{sample_id}.combined-cnv-calls.tsv"),
                  sample_id = all_samples)
  # Target PennCNV
  if TARGET == 'PennCNV':
    return expand(os.path.join(DATAPATH, "{sample_id}", "{sample_id}.penncnv-{chr}.tsv"),
                  sample_id = all_samples, chr = ('auto', 'chrx')) + \
      [os.path.join(DATAPATH, f"{sample_id}", f"{sample_id}.penncnv-chry.tsv") for
       sample_id, _, _, sex, _ in sample_data if sex.lower()[0] == 'm']
  # Target CBS
  if TARGET == 'CBS':
    return expand(os.path.join(DATAPATH, "{sample_id}", "{sample_id}.CBS.tsv"),
                  sample_id = all_samples)
  # Target CBS
  if TARGET == 'GADA':
    return expand(os.path.join(DATAPATH, "{sample_id}", "{sample_id}.GADA.tsv"),
                  sample_id = all_samples)
  # #Target: filtered-data
  # if TARGET == 'filtered-data':
  #   return  expand(os.path.join(DATAPATH, "{sample_id}", "{sample_id}.filtered-data.{filter}.tsv"),
  #                 sample_id = all_samples, filter = use_filter)

  #Target: unfiltered-data
  if TARGET == 'SNP-probe-data':
    return  expand(os.path.join(DATAPATH, "{sample_id}", "{sample_id}.processed-data.tsv"), sample_id = all_samples) + \
      expand(os.path.join(DATAPATH, "{sample_id}", "{sample_id}.unprocessed.vcf"), sample_id = all_samples)


rule all:
  input:
    get_target_files()
  run:
    if removetempconfig:
      os.remove(CONFIGFILE)


# Maybe-TODO / Alternative:
# - run in temp folder, then move all the individual gtc files (will need multi line bash script in that case, but that'd be fine)
rule run_gencall:
  input:
    bpm=config['static_data']['bpm_manifest_file'],
    egt=config['static_data']['egt_cluster_file'],
    idat_path = os.path.join(IDAT_INPUT, "{sentrix_name}")
  output:
    os.path.join(DATAPATH, "gtc", "{sentrix_name}", "_done")
  threads: get_tool_resource('GenCall', 'threads')
  resources:
    time=get_tool_resource('GenCall', 'runtime'),
    mem_mb=get_tool_resource('GenCall', 'memory'),
    partition=get_tool_resource('GenCall', 'partition')
  params:
    # TODO: Param chnages in this rule do not seem to trigger a rerun?
    options = get_tool_resource('GenCall', 'cmd-line-params'),
    outpath = os.path.join(DATAPATH, "gtc", "{sentrix_name}")
  log:
    err=os.path.join(LOGPATH, "GenCall", "{sentrix_name}", "error.log"),
    out=os.path.join(LOGPATH, "GenCall", "{sentrix_name}", "out.log")
  #TODO - containerize?
  # this wont ever be in conda, BUT there are docker containers (or even github repos) with this available
  # -> not clear if we should use a random docker for this (none by illumina) or make our own? Are we allowed even?
  shell:
    #TODO: check how imporant definition of the ICU version & LANG is here
    # CLR_ICU_VERSION_OVERRIDE="70.1" / $(uconv -V | sed 's/.* //g')
    'LANG="en_US.UTF-8" iaap-cli gencall "{input.bpm}" "{input.egt}" "{params.outpath}" --idat-folder "{input.idat_path}" --output-gtc {params.options} -t {threads} > {log.out} 2> {log.err} && [ $(ls {params.outpath}/*.gtc -l | wc -l) -ge 1 ] && touch {output} || exit 1'


# The iaap-cli will *always* generate filenames derived from the sentrix name & pos
def get_chip(wildcards, outtype = 'dir_path'):
  """Get the chip name from a sample_id
  Values for outtype: 'dirpath' | 'file'"""
  chip_name, chip_pos = [(n, p) for sid, n, p, _, _ in sample_data if sid == wildcards.sample_id][0]
  if outtype == 'dir_path':
    return os.path.join(DATAPATH, 'gtc', chip_name)
  elif outtype == 'file':
    return os.path.join(chip_name, chip_name + '_' + chip_pos + '.gtc')

#Note:
# the *.gtc output files from gencall will *always* match the idat_file names
# --> maybe better to switch from Chip_Name & Chip_Pos to Folder_Name & IDAT_Name ?!
rule relink_gencall:
  input:
    lambda wildcards: os.path.join(get_chip(wildcards), '_done')
  output:
    os.path.join(DATAPATH, "{sample_id}", "{sample_id}.gencall.gtc")
  params:
    gtc_link_path = lambda wildcards: os.path.join('..', 'gtc', get_chip(wildcards, outtype='file'))
  shell:
    'ln -s "{params.gtc_link_path}" "{output}"'
  

rule run_gtc2vcf_tsv:
  input:
    bpm=config['static_data']['bpm_manifest_file'],
    egt=config['static_data']['egt_cluster_file'],
    genome=config['static_data']['genome_fasta_file'],
    gtc = os.path.join(DATAPATH, "{sample_id}", "{sample_id}.gencall.gtc"),
  output:
    tsv = os.path.join(DATAPATH, "{sample_id}", "{sample_id}.processed-data.tsv"),
  threads: get_tool_resource('gtc2vcf', 'threads')
  resources:
    time=get_tool_resource('gtc2vcf', 'runtime'),
    mem_mb=get_tool_resource('gtc2vcf', 'memory'),
    partition=get_tool_resource('gtc2vcf', 'partition')
  params:
    options = get_tool_resource('gtc2vcf', 'cmd-line-params'),
    csv='--csv "{}"'.format(config['static_data']['csv_manifest_file']) if config['static_data']['csv_manifest_file'] else '',
  log:
    err=os.path.join(LOGPATH, "gtc2vcf", "{sample_id}", "tsv.error.log"),
    out=os.path.join(LOGPATH, "gtc2vcf", "{sample_id}", "tsv.out.log")
  shell:
    'bcftools plugin gtc2vcf {params.options} --no-version -O t --bpm "{input.bpm}" {params.csv} --egt "{input.egt}" --fasta-ref "{input.genome}" -o {output.tsv} {input.gtc} > {log.out} 2> {log.err}'

rule run_gtc2vcf_vcf:
  input:
    bpm=config['static_data']['bpm_manifest_file'],
    egt=config['static_data']['egt_cluster_file'],
    genome=config['static_data']['genome_fasta_file'],
    gtc = os.path.join(DATAPATH, "{sample_id}", "{sample_id}.gencall.gtc")
  output:
    vcf = os.path.join(DATAPATH, "{sample_id}", "{sample_id}.unprocessed.vcf"),
    metatxt = os.path.join(DATAPATH, "{sample_id}", "{sample_id}.stats.txt"),
  threads: get_tool_resource('gtc2vcf', 'threads')
  resources:
    time=get_tool_resource('gtc2vcf', 'runtime'),
    mem_mb=get_tool_resource('gtc2vcf', 'memory'),
    partition=get_tool_resource('gtc2vcf', 'partition')
  params:
    options = get_tool_resource('gtc2vcf', 'cmd-line-params'),
    csv='--csv "{}"'.format(config['static_data']['csv_manifest_file']) if config['static_data']['csv_manifest_file'] else '',
  log:
    err=os.path.join(LOGPATH, "gtc2vcf", "{sample_id}", "vcf.error.log"),
    out=os.path.join(LOGPATH, "gtc2vcf", "{sample_id}", "vcf.out.log"),
  shell: 'bcftools plugin gtc2vcf {params.options} -O v --bpm "{input.bpm}" {params.csv} --egt "{input.egt}" --fasta-ref "{input.genome}" --extra {output.metatxt} -o {output.vcf} {input.gtc} > {log.out} 2> {log.err}'

rule run_filter_tsv:
  input:
    tsv = os.path.join(DATAPATH, "{sample_id}", "{sample_id}.processed-data.tsv")
  output:
    tsv = os.path.join(DATAPATH, "{sample_id}", "{sample_id}.filtered-data-{filter}.tsv")
  threads: get_tool_resource('filter_tsv', 'threads')
  resources:
    time=get_tool_resource('filter_tsv', 'runtime'),
    mem_mb=get_tool_resource('filter_tsv', 'memory'),
    partition=get_tool_resource('filter_tsv', 'partition')
  log:
    err=os.path.join(LOGPATH, "filter_tsv", "{sample_id}", "{filter}.error.log"),
    out=os.path.join(LOGPATH, "filter_tsv", "{sample_id}", "{filter}.out.log")
  shell:
    'Rscript {SNAKEDIR}/scripts/filter_data.R -f {wildcards.filter} {input.tsv} {output.tsv} {CONFIGFILE} > {log.out} 2> {log.err}'
    
rule run_PennCNV:
  input:
    tsv = os.path.join(DATAPATH, "{sample_id}", "{{sample_id}}.filtered-data-{0}.tsv".format(get_tool_filter_settings('PennCNV'))),
    pfb = config['static_data']['pfb_file'],
    gcmodel = config['static_data']['GCmodel_file']
  output:
    tsv = os.path.join(DATAPATH, "{sample_id}", "{sample_id}.penncnv-{chr}.tsv"),
    err=os.path.join(LOGPATH,"PennCNV","{sample_id}","{chr}.error.log"),
  threads: get_tool_resource('PennCNV', 'threads')
  resources:
    time=get_tool_resource('PennCNV', 'runtime'),
    mem_mb=get_tool_resource('PennCNV', 'memory'),
    partition=get_tool_resource('PennCNV', 'partition')
  wildcard_constraints:
    chr='chrx|chry|auto'
  params:
    penncnv_path = '/home/user/PennCNV/detect_cnv.pl' if config['use_singularity'] else 'PennCNV_detect',
    filter = get_tool_filter_settings('PennCNV'),
    chrom = lambda wildcards: '' if wildcards.chr == 'auto' else '-'+wildcards.chr,
    #Male sex chromosomes can't have LOH, but PennCNV does not exclude it if run with -loh
    do_loh = lambda wildcards: '' if (wildcards.chr != 'auto' and get_ref_id(wildcards, True)[2] == 'm') else '-loh'
  log:
    err=os.path.join(LOGPATH, "PennCNV", "{sample_id}", "{chr}.error.log"),
    out=os.path.join(LOGPATH, "PennCNV", "{sample_id}", "{chr}.out.log")
  container:
    "docker://genomicslab/penncnv"
  shell:
    '{params.penncnv_path} -test {params.do_loh} {params.chrom} -confidence -hmm {SNAKEDIR}/PennCNV_overrides/hhall_loh.hmm -pfb {input.pfb} -gcmodel {input.gcmodel} {input.tsv} -out {output.tsv} > {log.out} 2> {log.err}'
    #'PennCNV_detect -test {params.do_loh} {params.chrom} -confidence -hmm {SNAKEDIR}/PennCNV_overrides/hhall_loh.hmm -pfb {input.pfb} -gcmodel {input.gcmodel} {input.tsv} -out {output.tsv} > {log.out} 2> {log.err}'


rule run_CBS:
  input:
    tsv = os.path.join(DATAPATH, "{sample_id}", "{{sample_id}}.filtered-data-{0}.tsv".format(get_tool_filter_settings('CBS')))
  output:
    tsv = os.path.join(DATAPATH, "{sample_id}", "{sample_id}.CBS.tsv")
  threads: get_tool_resource('CBS', 'threads')
  resources:
    time=get_tool_resource('CBS', 'runtime'),
    mem_mb=get_tool_resource('CBS', 'memory'),
    partition=get_tool_resource('CBS', 'partition')
  params:
    SDundo = config['settings']['CBS']['SDundo'],
    filter=get_tool_filter_settings('CBS'),
    settings=config['settings']['CBS']
  log:
    err=os.path.join(LOGPATH, "CBS", "{sample_id}", "error.log"),
    out=os.path.join(LOGPATH, "CBS", "{sample_id}", "out.log")
  shell:
    "Rscript {SNAKEDIR}/scripts/run_CBS_DNAcopy.R -s {params.SDundo} {input.tsv} {output.tsv} {CONFIGFILE} {SAMPLETABLE} > {log.out} 2> {log.err}"

rule run_GADA:
  input:
    tsv = os.path.join(DATAPATH, "{sample_id}", "{{sample_id}}.filtered-data-{0}.tsv".format(get_tool_filter_settings('GADA')))
  output:
    tsv = os.path.join(DATAPATH, "{sample_id}", "{sample_id}.GADA.tsv")
  threads: get_tool_resource('GADA', 'threads')
  resources:
    time=get_tool_resource('GADA', 'runtime'),
    mem_mb=get_tool_resource('GADA', 'memory'),
    partition=get_tool_resource('GADA', 'partition')
  params:
    filter = get_tool_filter_settings('GADA'),
    alpha = config['settings']['GADA']['aAlpha'],
    gadaT = config['settings']['GADA']['gada.T'],
    madT = config['settings']['GADA']['mad.T'],
    minLen = config['settings']['GADA']['MinSegLen']
  log:
    err=os.path.join(LOGPATH, "GADA", "{sample_id}", "err.log"),
    out=os.path.join(LOGPATH, "GADA", "{sample_id}", "out.log")
  shell:
    "Rscript {SNAKEDIR}/scripts/run_GADA.R -a {params.alpha} -g {params.gadaT} -m {params.madT} -l {params.minLen} {input.tsv} {output.tsv} > {log.out} 2> {log.err}"


def get_ref_id(wildcards, get_sex=False):
  sample_id = wildcards.sample_id
  sex, ref_id = [(s, r) for sid, _, _, s, r in sample_data if sid == sample_id][0]
  sex = sex[0].lower()

  if ref_id:
    try:
      #Assume existing match -> wrapper should have done a check
      ref_sex = [s for sid, _, _, s, _ in sample_data if sid == ref_id][0]
      ref_sex = ref_sex[0].lower()
    except IndexError:
      # Somehow no match
      raise SampletableReferenceError(f"Listed reference sample can not be found in sample-table: '{ref_id}'")
  else:
    ref_id = False
    ref_sex = False

  if get_sex:
    return sample_id, ref_id, sex, ref_sex
  else:
    return sample_id, ref_id

def get_preprocess_input(wildcards):
  sample_id, ref_id, sex, ref_sex = get_ref_id(wildcards, True)
  tools = config['settings']['CNV.calling.tools']
  files = []
  files += [os.path.join(DATAPATH, f"{sample_id}", f"{sample_id}.penncnv-auto.tsv"),
            os.path.join(DATAPATH, f"{sample_id}", f"{sample_id}.penncnv-chrx.tsv")] if 'PennCNV' in tools else []
  files += [os.path.join(DATAPATH, f"{sample_id}", f"{sample_id}.penncnv-chry.tsv")] if 'PennCNV' in tools and sex == 'm' else []
  files += [os.path.join(DATAPATH, f"{sample_id}", f"{sample_id}.CBS.tsv")] if 'CBS' in tools else []
  files += [os.path.join(DATAPATH, f"{sample_id}", f"{sample_id}.GADA.tsv")] if 'GADA' in tools else []

  files += [os.path.join(DATAPATH, f"{ref_id}", f"{ref_id}.combined-cnv-calls.tsv")] if ref_id else []

  return files

rule run_process_CNV_calls:
  input:
    get_preprocess_input
  output:
    os.path.join(DATAPATH, "{sample_id}", "{sample_id}.combined-cnv-calls.tsv")
  threads: get_tool_resource('CNV.process', 'threads')
  resources:
    time=get_tool_resource('CNV.process', 'runtime'),
    mem_mb=get_tool_resource('CNV.process', 'memory'),
    partition=get_tool_resource('CNV.process', 'partition')
  params:
    gada = '-g' if 'GADA' in config['settings']['CNV.calling.tools'] else '',
    penncnv = '-p' if 'PennCNV' in config['settings']['CNV.calling.tools'] else '',
    cbs = '-c' if 'CBS' in config['settings']['CNV.calling.tools'] else '',
    settings=config['settings']['postprocessing']
  log:
    err=os.path.join(LOGPATH, "CNV_process", "{sample_id}", "error.log"),
    out=os.path.join(LOGPATH, "CNV_process", "{sample_id}", "out.log")
  shell:
    "Rscript {SNAKEDIR}/scripts/process_CNV_calls.R {params.penncnv} {params.gada} {params.cbs} {DATAPATH} {wildcards.sample_id} {CONFIGFILE} {SAMPLETABLE} > {log.out} 2> {log.err}"


def get_report_sample_input(wildcards):
  sample_id, ref_id, sex, ref_sex = get_ref_id(wildcards, True)
  report_settings = config['reports'][wildcards.report]

  sample_files = expand(
          [os.path.join(DATAPATH, "{ids}", "{ids}.combined-cnv-calls.tsv"),
           os.path.join(DATAPATH, "{ids}", "{ids}.stats.txt"),
           os.path.join(LOGPATH, "PennCNV", "{ids}","{chrs}.error.log"),
           os.path.join(DATAPATH, "{ids}", "{ids}.processed-data.tsv"),
           os.path.join(DATAPATH, "{ids}", "{ids}.filtered-data-{filter}.tsv")],
          ids = (sample_id, ref_id) if ref_id else (sample_id,),
          chrs = ['auto', 'chrx'] + (['chry'] if sex == 'm' else []),
          filter = get_tool_filter_settings(f"report:{wildcards.report}:call.data.and.plots"),
  )

  incl_sections = config_extract(('include_sections', ), report_settings, config['reports']['__default__'])
  excl_sections = config_extract(('exclude_sections', ), report_settings, config['reports']['__default__'])
  do_snp_clustering = (incl_sections == '__all__' or 'SNP.dendrogram' in incl_sections) and not 'SNP.dendrogram' in excl_sections
  if do_snp_clustering:
    extra_sample_def = config_extract(('SNP_comparison', 'extra_samples'), report_settings, config['reports']['__default__'])
    ids = collect_SNP_cluster_ids(sample_id, extra_sample_def, sample_data_full)
    if not config_extract(('SNP_comparison', 'ignore_filter'), report_settings, config['reports']['__default__']):
      sample_files += expand([os.path.join(DATAPATH,"{ids}","{ids}.filtered-data-{filter}.tsv")],
        ids=ids, filter = get_tool_filter_settings(f"report:{wildcards.report}:SNP_comparison"))
    else:
      sample_files += expand([os.path.join(DATAPATH, "{ids}", "{ids}.processed-data.tsv")], ids = ids)

  return sample_files

rule knit_report:
  input:
    get_report_sample_input
  output:
    report=os.path.join(DATAPATH, "{sample_id}", "{sample_id}.{report}.{ext}") ,
    plots=directory(os.path.join(DATAPATH, "{sample_id}", "{report}-{ext}_images")),
    # -> this extra output is optional now!
    #summary=os.path.join(DATAPATH, "{sample_id}", "{sample_id}.summary-check.tsv")
  wildcard_constraints:
    ext="(pdf|html)"
  threads: get_tool_resource('knitr', 'threads')
  resources:
    time=get_tool_resource('knitr', 'runtime'),
    mem_mb=get_tool_resource('knitr', 'memory'),
    partition=get_tool_resource('knitr', 'partition')
  log:
    err=os.path.join(LOGPATH, "report", "{sample_id}", "{report}-{ext}.error.log"),
    out=os.path.join(LOGPATH, "report", "{sample_id}", "{report}-{ext}.out.log")
  shell:
    "Rscript {SNAKEDIR}/scripts/knit_report.R {wildcards.sample_id} {wildcards.report} {CONFIGFILE} > {log.out} 2> {log.err}"


rule run_make_cnv_vcf:
  input:
    os.path.join(DATAPATH,"{sample_id}","{sample_id}.combined-cnv-calls.tsv"),
    os.path.join(DATAPATH,"{sample_id}","{sample_id}.unprocessed.vcf"),
  output:
    get_cnv_vcf_output(config['settings']['make_cnv_vcf']['mode'])
  threads: get_tool_resource('make_cnv_vcf', 'threads')
  resources:
    time=get_tool_resource('make_cnv_vcf', 'runtime'),
    mem_mb=get_tool_resource('make_cnv_vcf', 'memory'),
    partition=get_tool_resource('make_cnv_vcf', 'partition')
  params:
    include_states=' '.join(config['settings']['make_cnv_vcf']['include_states']),
    mode = config['settings']['make_cnv_vcf']['mode']
  log:
    err=os.path.join(LOGPATH,"make_cnv_vcf","{sample_id}","error.log"),
    out=os.path.join(LOGPATH,"make_cnv_vcf","{sample_id}","out.log")
  shell:
    "Rscript {SNAKEDIR}/scripts/make_cnv_vcf.R {DATAPATH} {wildcards.sample_id} {CONFIGFILE} -m {params.mode} -i {params.include_states} > {log.out} 2> {log.err}"
