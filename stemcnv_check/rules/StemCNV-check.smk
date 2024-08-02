# -*- coding: utf-8 -*-

import importlib.resources
import os
from pathlib import Path
import tempfile
import ruamel.yaml as ruamel_yaml
from stemcnv_check import STEM_CNV_CHECK
from stemcnv_check.helpers import read_sample_table, config_extract, collect_SNP_cluster_ids
from stemcnv_check.exceptions import SampleConstraintError, ConfigValueError

SNAKEDIR = str(importlib.resources.files(STEM_CNV_CHECK))

# Configuration ================================================================
if not config:
  # This should not happen unless the snakefile is invoked directly from cmd-line
  # try to load a baseline config file
  CONFIGFILE = "config.yaml"
  configfile: CONFIGFILE
  removetempconfig = False
else:
  # This a manual workaround around for passing the config as a snakemake object to R & python scripts
  # Save config to tempfile & use that instead, this ensures the combined default + user configs are used
  # It also archives all options passed by command line and allows them all to be displayed in report
  f, CONFIGFILE = tempfile.mkstemp(suffix = '.yaml', text=True)
  os.close(f)
  removetempconfig = True

yaml = ruamel_yaml.YAML(typ='safe')
config['snakedir'] = SNAKEDIR
with open(CONFIGFILE, 'w') as yamlout:
  yaml.dump(config, yamlout)


SAMPLETABLE = config['sample_table'] if 'sample_table' in config else 'sample_table.tsv' # Defined by wrapper
BASEPATH = config['basedir'] if 'basedir' in config else os.getcwd() #Defined by wrapper
DATAPATH = config['data_path'] #if os.path.isabs(config['data_path']) else os.path.join(BASEPATH, config['data_path'])
LOGPATH = config['log_path'] #if os.path.isabs(config['log_path']) else os.path.join(BASEPATH, config['log_path'])
TARGET = config['target'] if 'target' in config else 'report' #Defined by wrapper
IDAT_INPUT = config['raw_data_folder']

wildcard_constraints:
  sample_id=config['wildcard_constraints']['sample_id'],
  #sentrix_name=config['wildcard_constraints']['sentrix_name'],
  #sentrix_pos=config['wildcard_constraints']['sentrix_pos'],

#Never submit these to cluster
localrules:
  relink_gencall,
  all

sample_data = read_sample_table(SAMPLETABLE)
sample_data_full = read_sample_table(SAMPLETABLE, with_opt=True)

include: 'common.smk'
include: 'illumina_raw_processing.smk'
include: "SNP_processing.smk"

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
  #Target options: ('report', 'cnv-vcf', 'combined-cnv-calls', 'PennCNV', 'CBS', 'SNP-probe-data'),
  all_samples = [sample_id for sample_id, _, _, _, _ in sample_data]

  # complete
  if TARGET == 'complete':
    return expand(os.path.join(DATAPATH, "{sample_id}", "{sample_id}.{report_filetype}"),
                  sample_id = all_samples,
                  report_filetype = [rep+'.'+config['reports'][rep]['file_type'] for rep in config['reports'].keys() if rep != '__default__']) + \
           expand(get_cnv_vcf_output(config['settings']['make_cnv_vcf']['mode']),
                  sample_id = all_samples)
  # expand(os.path.join(DATAPATH,"{sample_id}","{sample_id}.summary-check.tsv"),
  #        sample_id = all_samples)

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

  # Target Processed-calls
  if TARGET == 'combined-cnv-calls':
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
  # #Target: filtered-data
  # if TARGET == 'filtered-data':
  #   return  expand(os.path.join(DATAPATH, "{sample_id}", "{sample_id}.filtered-data.{filter}.tsv"),
  #                 sample_id = all_samples, filter = use_filter)

  #Target: unfiltered-data
  if TARGET == 'SNP-probe-data':
    return expand(os.path.join(DATAPATH, "{sample_id}", "{sample_id}.processed-data.tsv"), sample_id = all_samples) + \
      expand(os.path.join(DATAPATH, "{sample_id}", "{sample_id}.unprocessed.vcf"), sample_id = all_samples)

rule all:
  input:
    get_target_files()
  run:
    if removetempconfig:
      os.remove(CONFIGFILE)


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
    runtime=get_tool_resource('gtc2vcf', 'runtime'),
    mem_mb=get_tool_resource('gtc2vcf', 'memory'),
    partition=get_tool_resource('gtc2vcf', 'partition')
  params:
    options = get_tool_resource('gtc2vcf', 'cmd-line-params'),
    csv='--csv "{}"'.format(config['static_data']['csv_manifest_file']) if config['static_data']['csv_manifest_file'] else '',
  log:
    err=os.path.join(LOGPATH, "gtc2vcf", "{sample_id}", "tsv.error.log"),
    out=os.path.join(LOGPATH, "gtc2vcf", "{sample_id}", "tsv.out.log")
  conda:
    importlib.resources.files(STEM_CNV_CHECK).joinpath("envs","gtc2vcf.yaml")
  shell:
    'bcftools plugin gtc2vcf {params.options} --no-version -O t --bpm "{input.bpm}" {params.csv} --egt "{input.egt}" --fasta-ref "{input.genome}" -o {output.tsv} {input.gtc} > {log.out} 2> {log.err}'


# rule run_filter_tsv:
#   input:
#     tsv = os.path.join(DATAPATH, "{sample_id}", "{sample_id}.processed-data.tsv")
#   output:
#     tsv = os.path.join(DATAPATH, "{sample_id}", "{sample_id}.filtered-data-{filter}.tsv")
#   threads: get_tool_resource('filter_tsv', 'threads')
#   resources:
#     runtime=get_tool_resource('filter_tsv', 'runtime'),
#     mem_mb=get_tool_resource('filter_tsv', 'memory'),
#     partition=get_tool_resource('filter_tsv', 'partition')
#   log:
#     err=os.path.join(LOGPATH, "filter_tsv", "{sample_id}", "{filter}.error.log"),
#     out=os.path.join(LOGPATH, "filter_tsv", "{sample_id}", "{filter}.out.log")
#   conda:
#     importlib.resources.files(STEM_CNV_CHECK).joinpath("envs","general-R.yaml")
#   shell:
#     'Rscript {SNAKEDIR}/scripts/filter_data.R -f {wildcards.filter} {input.tsv} {output.tsv} {CONFIGFILE} > {log.out} 2> {log.err}'
# 
# 
# rule run_PennCNV:
#   input:
#     tsv=os.path.join(DATAPATH,"{sample_id}","{{sample_id}}.filtered-data-{0}.tsv".format(get_tool_filter_settings('PennCNV'))),
#     pfb=config['static_data']['penncnv_pfb_file'],
#     gcmodel=config['static_data']['penncnv_GCmodel_file']
#   output:
#     tsv=os.path.join(DATAPATH,"{sample_id}","{sample_id}.penncnv-{chr}.tsv"),
#     err=os.path.join(LOGPATH,"PennCNV","{sample_id}","{chr}.error.log")
#   threads: get_tool_resource('PennCNV','threads')
#   resources:
#     runtime=get_tool_resource('PennCNV','runtime'),
#     mem_mb=get_tool_resource('PennCNV','memory'),
#     partition=get_tool_resource('PennCNV','partition')
#   wildcard_constraints:
#     chr='chrx|chry|auto'
#   params:
#     filter=get_tool_filter_settings('PennCNV'),
#     chrom=lambda wildcards: '' if wildcards.chr == 'auto' else '-' + wildcards.chr,
#     #Male sex chromosomes can't have LOH, but PennCNV does not exclude it if run with -loh
#     do_loh=lambda wildcards: '' if (wildcards.chr != 'auto' and get_ref_id(wildcards,True)[2] == 'm') else '-loh',
#     snakedir = fix_container_path(SNAKEDIR, 'snakedir'),
#     tsvout=lambda wildcards: fix_container_path(os.path.join(DATAPATH, wildcards.sample_id, wildcards.sample_id+".penncnv-"+wildcards.chr+".tsv"),'data'),
#     tsvin=lambda wildcards: fix_container_path(os.path.join(DATAPATH, wildcards.sample_id, wildcards.sample_id+".filtered-data-{}.tsv".format(get_tool_filter_settings('PennCNV'))),'data'),
#     pfb=fix_container_path(config['static_data']['penncnv_pfb_file'],'static'),
#     gcmodel=fix_container_path(config['static_data']['penncnv_GCmodel_file'],'static'),
#     logerr=lambda wildcards:fix_container_path(os.path.join(LOGPATH,"PennCNV", wildcards.sample_id, wildcards.chr+".error.log"),'logs'),
#     logout=lambda wildcards:fix_container_path(os.path.join(LOGPATH,"PennCNV", wildcards.sample_id, wildcards.chr+".out.log"),'logs')
#   log:
#     err=os.path.join(LOGPATH,"PennCNV","{sample_id}","{chr}.error.log"),
#     out=os.path.join(LOGPATH,"PennCNV","{sample_id}","{chr}.out.log")
#   container:
#     "docker://genomicslab/penncnv"
#   shell:
#     '/home/user/PennCNV/detect_cnv.pl -test {params.do_loh} {params.chrom} -confidence -hmm {params.snakedir}/supplemental-files/hhall_loh.hmm -pfb {params.pfb} -gcmodel {params.gcmodel} {params.tsvin} -out {params.tsvout} > {params.logout} 2> {params.logerr}'


rule run_CBS:
  input:
    tsv = os.path.join(DATAPATH, "{sample_id}", "{{sample_id}}.filtered-data-{0}.tsv".format(get_tool_filter_settings('CBS')))
  output:
    tsv = os.path.join(DATAPATH, "{sample_id}", "{sample_id}.CBS.tsv")
  threads: get_tool_resource('CBS', 'threads')
  resources:
    runtime=get_tool_resource('CBS', 'runtime'),
    mem_mb=get_tool_resource('CBS', 'memory'),
    partition=get_tool_resource('CBS', 'partition')
  params:
    SDundo = config['settings']['CBS']['SDundo'],
    filter=get_tool_filter_settings('CBS'),
    settings=config['settings']['CBS']
  log:
    err=os.path.join(LOGPATH, "CBS", "{sample_id}", "error.log"),
    out=os.path.join(LOGPATH, "CBS", "{sample_id}", "out.log")
  conda:
    importlib.resources.files(STEM_CNV_CHECK).joinpath("envs", "general-R.yaml")
  shell:
    "Rscript {SNAKEDIR}/scripts/run_CBS_DNAcopy.R -s {params.SDundo} {input.tsv} {output.tsv} {CONFIGFILE} {SAMPLETABLE} > {log.out} 2> {log.err}"


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
      raise SampleConstraintError(f"Listed reference sample can not be found in sample-table: '{ref_id}'")
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

  filter = get_tool_filter_settings(f"settings:CNV_processing:call_processing")
  files += [os.path.join(DATAPATH, f"{sample_id}", f"{sample_id}.filtered-data-{filter}.tsv")]

  files += [os.path.join(DATAPATH, f"{ref_id}", f"{ref_id}.combined-cnv-calls.tsv")] if ref_id else []

  return files

rule run_process_CNV_calls:
  input:
    get_preprocess_input
  output:
    os.path.join(DATAPATH, "{sample_id}", "{sample_id}.combined-cnv-calls.tsv")
  threads: get_tool_resource('CNV.process', 'threads')
  resources:
    runtime=get_tool_resource('CNV.process', 'runtime'),
    mem_mb=get_tool_resource('CNV.process', 'memory'),
    partition=get_tool_resource('CNV.process', 'partition')
  params:
    penncnv = '-p' if 'PennCNV' in config['settings']['CNV.calling.tools'] else '',
    cbs = '-c' if 'CBS' in config['settings']['CNV.calling.tools'] else '',
    settings=config['settings']['CNV_processing']
  log:
    err=os.path.join(LOGPATH, "CNV_process", "{sample_id}", "error.log"),
    out=os.path.join(LOGPATH, "CNV_process", "{sample_id}", "out.log")
  conda:
    importlib.resources.files(STEM_CNV_CHECK).joinpath("envs","general-R.yaml")
  shell:
    "Rscript {SNAKEDIR}/scripts/process_CNV_calls.R {params.penncnv} {params.cbs} {DATAPATH} {wildcards.sample_id} {CONFIGFILE} {SAMPLETABLE} > {log.out} 2> {log.err}"


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

  if wildcards.ext == 'pdf':
    sample_files += [os.path.join(LOGPATH, "report", "_latex_installation_check")]

  return sample_files

rule check_latex_installation:
  output:
    os.path.join(LOGPATH, "report", "_latex_installation_check")
  conda:
    importlib.resources.files(STEM_CNV_CHECK).joinpath("envs","general-R.yaml")
  shell:
    """
Rscript - << 'EOF'
suppressMessages(library(tinytex))

latex_path <- tinytex_root()

if (latex_path == "") {{
  message("No LaTeX installation found. Installing TinyTex.")
  install_tinytex()
}}

system("touch {output}")

EOF
"""

# def get_report_config(wildcards):
#   return config['reports'][wildcards.report]

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
  # params:
  #   config = get_report_config
  resources:
    runtime=get_tool_resource('knitr', 'runtime'),
    mem_mb=get_tool_resource('knitr', 'memory'),
    partition=get_tool_resource('knitr', 'partition')
  log:
    err=os.path.join(LOGPATH, "report", "{sample_id}", "{report}-{ext}.error.log"),
    out=os.path.join(LOGPATH, "report", "{sample_id}", "{report}-{ext}.out.log")
  conda:
    importlib.resources.files(STEM_CNV_CHECK).joinpath("envs","general-R.yaml")
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
    runtime=get_tool_resource('make_cnv_vcf', 'runtime'),
    mem_mb=get_tool_resource('make_cnv_vcf','memory'),
    partition=get_tool_resource('make_cnv_vcf', 'partition')
  params:
    include_states=' '.join(config['settings']['make_cnv_vcf']['include_states']),
    mode = config['settings']['make_cnv_vcf']['mode']
  log:
    err=os.path.join(LOGPATH,"make_cnv_vcf","{sample_id}","error.log"),
    out=os.path.join(LOGPATH,"make_cnv_vcf","{sample_id}","out.log")
  conda:
    importlib.resources.files(STEM_CNV_CHECK).joinpath("envs","general-R.yaml")
  shell:
    "Rscript {SNAKEDIR}/scripts/make_cnv_vcf.R {DATAPATH} {wildcards.sample_id} {CONFIGFILE} -m {params.mode} -i {params.include_states} > {log.out} 2> {log.err}"
