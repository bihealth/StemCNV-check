# -*- coding: utf-8 -*-

import os
from scripts.py_helpers import read_sample_table
# Configuration ================================================================

#ONLY load if somehow manually started
if not config:
	CONFIGFILE = "config.yaml"
	configfile: CONFIGFILE
else:
	CONFIGFILE = config['configfile']

SAMPLETABLE = config['sample_table'] if 'sample_table' in config else 'sample_table.txt' # Defined by wrapper
SNAKEDIR = config['snakedir'] if 'snakedir' in config else os.getcwd() #Defined by wrapper 
BASEPATH = config['basedir'] if 'basedir' in config else os.getcwd() #Defined by wrapper 
IDAT_INPUT = config['raw_data_folder']

#TODO enforcin chip_no & pos could be problematic for people?
# -> think of a more flexible way to define files ?
#TODO load this from config
wildcard_constraints: 
	chr=config['wildcard_constraints']['chr'],
	filter=config['wildcard_constraints']['filter'],
	sample_id=config['wildcard_constraints']['sample_id'],
	chip_name=config['wildcard_constraints']['chip_name'],
	
#Never submit these to cluster
localrules:
	relink_gencall,

sample_data = read_sample_table(SAMPLETABLE)

# Rules ========================================================================

rule all:
  input: 
    # Final / Report
    [os.path.join(BASEPATH, "data", f"{sample_id}", f"{sample_id}.CNV-report.html") for
      _, sample_id, _, _ in sample_data.values()] 
    # # Target PennCNV
    # [os.path.join(BASEPATH, "data", f"{sample_id}", f"{sample_id}.penncnv-autosomes.tsv") for
    #   _, sample_id, _, _ in sample_data.values()] + 
    # [os.path.join(BASEPATH, "data", f"{sample_id}", f"{sample_id}.penncnv-chrx.tsv") for
    #   _, sample_id, _, _ in sample_data.values()] +  
    # [os.path.join(BASEPATH, "data", f"{sample_id}", f"{sample_id}.penncnv-chry.tsv") for
    #   _, sample_id, _, _ in sample_data.values() if sex == 'Male'] + 
    # # Target CBS
    # [os.path.join(BASEPATH, "data", f"{sample_id}", f"{sample_id}.CBS.tsv") for
    #   chip_no, chip_pos, _, _ in sample_data.values()] +
    # #Target: processed data
    # [os.path.join(BASEPATH, "data", f"{sample_id}", f"{sample_id}.processed-data.tsv") for
    #   chip_no, chip_pos, _, _ in sample_data.values()] +
    # [os.path.join(BASEPATH, "data", f"{sample_id}", f"{sample_id}.unprocessed.vcf") for
    #   chip_no, chip_pos, _, _ in sample_data.values()]
    

# TODO: make he gtc files intermediate? tehy are only needed once really
# -> would allow an easier structure for output data
# Optional(/Alternative?) TODO:
# - make temp folder that link to specific idat files only so that gencall is only run a one a sample at a time
rule run_gencall:
  input:
    bpm=config['static_data']['bpm_manifest_file'],
    egt=config['static_data']['egt_cluster_file'],
    idat_path = os.path.join(IDAT_INPUT, "{chip_name}")
  output:
    os.path.join(BASEPATH, "data", "gtc", "{chip_name}", "_done")
  threads: config['tools']['GenCall']['threads']
  resources:
    time=config['tools']['GenCall']['runtime'],
    mem_mb=config['tools']['GenCall']['memory'],
    partition='medium'
  params:
    #TODO: check which options make sense here
    options = config['tools']['GenCall']['cmd-line-params'],
    outpath = os.path.join(BASEPATH, "data", "gtc", "{chip_name}")
  log:
    err=os.path.join(BASEPATH, "logs", "GenCall", "{chip_name}", "error.log"),
    out=os.path.join(BASEPATH, "logs", "GenCall", "{chip_name}", "out.log")
  shell: 
    #TODO: check how imporant definition of the ICU version & LANG is here
    # CLR_ICU_VERSION_OVERRIDE="70.1" / $(uconv -V | sed 's/.* //g')
    'LANG="en_US.UTF-8" iaap-cli gencall {input.bpm} {input.egt} {params.outpath} --idat-folder {input.idat_path} --output-gtc {params.options} -t {threads} > {log.out} 2> {log.err} && [ $(ls {params.outpath}/*.gtc -l | wc -l) -ge 1 ] && touch {output} || exit 1'



def get_gtc_path(wildcards):
	"Get the chip name from  a sample_id"
	chip_name = [n for _,(n, i, _, _) in sample_data.items() if i == wildcards.sample_id][0]
	return  os.path.join(BASEPATH, 'data', 'gtc', chip_name)
	
rule relink_gencall:
  input:
    lambda wildcards: os.path.join(get_gtc_path(wildcards), '_done')
  output:
    os.path.join(BASEPATH, "data", "{sample_id}", "{sample_id}.gencall.gtc")
  params:
    gtc_file = lambda wildcards: os.path.join(get_gtc_path(wildcards), f"{wildcards.sample_id}.gtc")
  shell:
    "ln -s {params.gtc_file} {output}"
  

rule run_gtc2vcf_tsv:
  input:
    bpm=config['static_data']['bpm_manifest_file'],
    egt=config['static_data']['egt_cluster_file'],
    genome=config['static_data']['genome_fasta_file'],
    gtc = os.path.join(BASEPATH, "data", "{sample_id}", "{sample_id}.gencall.gtc"),
  output:
    tsv = os.path.join(BASEPATH, "data", "{sample_id}", "{sample_id}.processed-data.tsv"),
  threads: 1
  resources:
    time=config['tools']['gtc2vcf']['runtime'],
    mem_mb=config['tools']['gtc2vcf']['memory'],
    partition='medium'
  params:
    options = config['tools']['gtc2vcf']['cmd-line-params'],
    csv='--csv {}'.format(config['static_data']['csv_manifest_file']) if config['static_data']['csv_manifest_file'] else '',
  log:
    err=os.path.join(BASEPATH, "logs", "gtc2vcf", "{sample_id}", "tsv.error.log"),
    out=os.path.join(BASEPATH, "logs", "gtc2vcf", "{sample_id}", "tsv.out.log")
  shell:
    'bcftools plugin gtc2vcf {params.options} --no-version -O t --bpm {input.bpm} {params.csv} --egt {input.egt} --fasta-ref {input.genome} -o {output.tsv} {input.gtc} > {log.out} 2> {log.err}'

rule run_gtc2vcf_vcf:
  input:
    bpm=config['static_data']['bpm_manifest_file'],
    egt=config['static_data']['egt_cluster_file'],
    genome=config['static_data']['genome_fasta_file'],
    gtc = os.path.join(BASEPATH, "data", "{sample_id}", "{sample_id}.gencall.gtc")
  output:
    vcf = os.path.join(BASEPATH, "data", "{sample_id}", "{sample_id}.unprocessed.vcf"),
    metatxt = os.path.join(BASEPATH, "data", "{sample_id}", "{sample_id}.stats.txt"),
  threads: 1
  resources:
    time=config['tools']['gtc2vcf']['runtime'],
    mem_mb=config['tools']['gtc2vcf']['memory'],
    partition='medium',
  params:
    options = config['tools']['gtc2vcf']['cmd-line-params'],
    csv='--csv {}'.format(config['static_data']['csv_manifest_file']) if config['static_data']['csv_manifest_file'] else '',
  log:
    err=os.path.join(BASEPATH, "logs", "gtc2vcf", "{sample_id}", "vcf.error.log"),
    out=os.path.join(BASEPATH, "logs", "gtc2vcf", "{sample_id}", "vcf.out.log"),
  shell:
    'bcftools plugin gtc2vcf {params.options} --no-version -O v --bpm {input.bpm} {params.csv} --egt {input.egt} --fasta-ref {input.genome} --extra {output.metatxt} -o {output.vcf} {input.gtc} > {log.out} 2> {log.err}'
 
rule filter_tsv:
  input:
    tsv = os.path.join(BASEPATH, "data", "{sample_id}", "{sample_id}.processed-data.tsv")
  output:
    tsv = os.path.join(BASEPATH, "data", "{sample_id}", "{sample_id}.filtered-data.{filter}.tsv")
  threads: 1
  log:
    err=os.path.join(BASEPATH, "logs", "filter_tsv", "{sample_id}", "{filter}.error.log"),
    out=os.path.join(BASEPATH, "logs", "filter_tsv", "{sample_id}", "{filter}.out.log")
  shell:
    'Rscript {SNAKEDIR}/scripts/filter_data.R -f {wildcards.filter} {input.tsv} {output.tsv} > {log.out} 2> {log.err}'
    
rule make_PennCNV_sexfile:
  output:
    filename=temp(os.path.join(BASEPATH, "penncnv-sexfile.txt"))
  # a 2-column file containing filename and sex (male/female) for
  # sex chromosome calling with -chrx argument. The first
  # tab-delimited column should be the input signal file name, while
  # the second tab-delimited column should be male or female.
  # Alternatively, abbreviations including m (male), f (female), 1
  # (male) or 2 (female) are also fine.
  params:
    filter = config['settings']['filter']['use-filterset']
  run:
    with open(output.filename, 'w') as f:
      for _, sample_id, sex, _ in sample_data.values():
        inputfile = os.path.join(BASEPATH, "data", f"{sample_id}", f"{sample_id}.filtered-data.{params.filter}.tsv")
        #ensure its consistently 'm'/'f'
        sex = sex.lower()[0]
        f.write(f"{inputfile}\t{sex}\n")

rule run_PennCNV_auto:
  input:
    tsv = os.path.join(BASEPATH, "data", "{sample_id}", "{sample_id}.filtered-data.{filter}.tsv"),
    pfb = config['static_data']['pfb_file'],
    gcmodel = config['static_data']['GCmodel_file']
  output:
    tsv = os.path.join(BASEPATH, "data", "{sample_id}", "{sample_id}.penncnv-autosomes.{filter}.tsv")
  threads: 1
  # resources:
  #   time=...
  #   memory=...
  #   partition=...
  log:
    err=os.path.join(BASEPATH, "logs", "PennCNV", "{sample_id}", "auto.{filter}.error.log"),
    out=os.path.join(BASEPATH, "logs", "PennCNV", "{sample_id}", "auto.{filter}.out.log")
  shell:
    'PennCNV_detect -test -loh -confidence -hmm {SNAKEDIR}/PennCNV_overrides/hhall_loh.hmm -pfb {input.pfb} -gcmodel {input.gcmodel} {input.tsv} -out {output.tsv} > {log.out} 2> {log.err}'


rule run_PennCNV_sex:
  input:
    tsv = os.path.join(BASEPATH, "data", "{sample_id}", "{sample_id}.filtered-data.{filter}.tsv"),
    sexfile = ancient(os.path.join(BASEPATH, "penncnv-sexfile.txt")),
    pfb = config['static_data']['pfb_file'],
    gcmodel = config['static_data']['GCmodel_file']
  output:
    tsv = os.path.join(BASEPATH, "data", "{sample_id}", "{sample_id}.penncnv-{chr}.{filter}.tsv")
  threads: 1
  # resources:
  #   time=...
  #   memory=...
  #   partition=...
  log:
    err=os.path.join(BASEPATH, "logs", "PennCNV", "{sample_id}", "{chr}.{filter}.error.log"),
    out=os.path.join(BASEPATH, "logs", "PennCNV", "{sample_id}", "{chr}.{filter}.out.log")
  shell:
    'PennCNV_detect -test -loh -confidence -hmm {SNAKEDIR}/PennCNV_overrides/hhall_loh.hmm -pfb {input.pfb} -gcmodel {input.gcmodel} -{wildcards.chr} -sex {input.sexfile} {input.tsv} -out {output.tsv} > {log.out} 2> {log.err}'

#? maybe ??
# make a rule to merge PennCNV files?
# -> that rule would also need to check sex though

rule run_CBS:
  input:
    tsv = os.path.join(BASEPATH, "data", "{sample_id}", "{sample_id}.filtered-data.{filter}.tsv")
  output:
    tsv = os.path.join(BASEPATH, "data", "{sample_id}", "{sample_id}.CBS.{filter}.tsv")
  threads: 1
  params:
    SDundo = config['settings']['CBS']['SDundo']
  log:
    err=os.path.join(BASEPATH, "logs", "CBS", "{sample_id}", "{filter}.error.log"),
    out=os.path.join(BASEPATH, "logs", "CBS", "{sample_id}", "{filter}.out.log")
  shell:
    "Rscript {SNAKEDIR}/scripts/run_CBS_DNAcopy.R -s {params.SDundo} {input.tsv} {output.tsv} > {log.out} 2> {log.err}"

# rule run_GADA:
#   another CNV algorithm; also mosaicism?
# 
# rule process_CNV_calls:
#   -> currently done in report; move out ?
   
def get_report_sample_input(wildcards):
  sex, ref_name = [(s, r) for name,(c, i, s, r) in sample_data.items() if i == wildcards.sample_id][0]
  sample_id = wildcards.sample_id
  if ref_name and ref_name in sample_data:
    _, ref_id, ref_sex, _ = sample_data[ref_name]
  elif ref_name:
    # listed reference does not exist in sampletable
    # TODO make a nicer exception class
    raise Exception(f"Listed reference sample can not be found in sample-table: '{ref_name}'")
  else:
    ref_id = False
    
  files = expand(
          [os.path.join(BASEPATH, "data", "{ids}", "{ids}.penncnv-autosomes.{filter}.tsv"),
           os.path.join(BASEPATH, "data", "{ids}", "{ids}.penncnv-chrx.{filter}.tsv"),
           os.path.join(BASEPATH, "data", "{ids}", "{ids}.CBS.{filter}.tsv"),
           os.path.join(BASEPATH, "data", "{ids}", "{ids}.stats.txt"),
           os.path.join(BASEPATH, "data", "{ids}", "{ids}.processed-data.tsv"),],
           ids = (sample_id, ref_id) if ref_id else (sample_id,),
           filter = config['settings']['filter']['use-filterset']
          )
  filter = config['settings']['filter']['use-filterset']
  files = files + [os.path.join(BASEPATH, "data", f"{sample_id}", f"{sample_id}.penncnv-chry.{filter}.tsv")] if sex[0].lower() == 'm' else files 
  files = files + [os.path.join(BASEPATH, "data", f"{ref_id}", f"{ref_id}.penncnv-chry.{filter}.tsv")] if (ref_id and ref_sex[0].lower() == 'm') else files

  return files


rule knit_report:
  input:
    get_report_sample_input
  output:
    html=os.path.join(BASEPATH, "data", "{sample_id}", "{sample_id}.CNV-report.html"),
    plots=directory(os.path.join(BASEPATH, "data", "{sample_id}", "report_images"))
  resources:
    time=config['tools']['knitr']['runtime'],
    mem_mb=config['tools']['knitr']['memory'],
    partition='medium'
  log:
    err=os.path.join(BASEPATH, "logs", "report", "{sample_id}", "error.log"),
    out=os.path.join(BASEPATH, "logs", "report", "{sample_id}", "out.log")
  shell:
    "Rscript {SNAKEDIR}/scripts/knit_report.R {wildcards.sample_id} {BASEPATH} {SNAKEDIR} {SAMPLETABLE} {CONFIGFILE} > {log.out} 2> {log.err}"
