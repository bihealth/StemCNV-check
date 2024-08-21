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
include: "penncnv.smk"

include: "report_generation.smk"

# Rule function,s move to common? ========================================================================

def get_cnv_vcf_output(mode):
    addition = config['settings']['make_cnv_vcf']['name_addition']
    if mode == 'combined-calls':
        return os.path.join(DATAPATH, "{sample_id}", f"{{sample_id}}.combined-cnv-calls{addition}.vcf")
    elif mode == 'split-tools':
        return [os.path.join(DATAPATH, "{sample_id}", f"{{sample_id}}.{tool}-cnv-calls{addition}.vcf")
                for tool in config['settings']['CNV.calling.tools']]
    else:
        raise ConfigValueError('Value not allowed for settings$make.cnv.vcf$mode: "{}"'.format(config['settings']['make_cnv_vcf']['mode']))

def get_target_files(target = TARGET):
    #Target options: ('report', 'combined-cnv-calls', 'PennCNV', 'CBS', 'SNP-data'),
    all_samples = [sample_id for sample_id, _, _, _, _ in sample_data]

    # complete
    if target == 'complete':
        return get_target_files('report') + get_target_files('combined-cnv-calls')

    # Report
    if target == 'report':
        return expand(
            os.path.join(DATAPATH, "{sample_id}", "{sample_id}.{report_filetype}"),
            sample_id = all_samples,
            report_filetype = [rep+'.'+config['reports'][rep]['file_type'] for rep in config['reports'].keys() 
                               if rep != '__default__']) #+ \
                     # expand(os.path.join(DATAPATH,"{sample_id}","{sample_id}.summary-check.tsv"),
                     #                sample_id = all_samples)
    # Target Processed-calls
    if target == 'combined-cnv-calls':
        return expand([
                os.path.join(DATAPATH, "{sample_id}", "{sample_id}.combined-cnv-calls.tsv"),
                os.path.join(DATAPATH, "{sample_id}", "{sample_id}.combined-cnv-calls.vcf.gz"),
            ],
            sample_id = all_samples
        )
    # Target PennCNV
    if target == 'PennCNV':
        return expand(
            os.path.join(DATAPATH, "{sample_id}", "{sample_id}.CNV_calls.penncnv.vcf.gz"),
            sample_id = all_samples,
        )
    # Target CBS
    if target == 'CBS':
        return expand(
            os.path.join(DATAPATH, "{sample_id}", "{sample_id}.CNV_calls.CBS.vcf.gz"),
            sample_id = all_samples
        )
    # Target SNP-data
    if target == 'SNP-data':
        #TODO: update this
        return expand(
            os.path.join(DATAPATH,"{sample_id}","{sample_id}.processed-SNP-data.{filter}-filter.vcf"),
            # os.path.join(DATAPATH, "{sample_id}", "{sample_id}.annotated-SNP-data.{filter}-filter.vcf.gz"),
            sample_id = all_samples,
            filter = config['settings']['default-filter-set']
        ) 
    
# Rules ========================================================================

rule all:
    input:
        get_target_files()
    run:
        if removetempconfig:
            os.remove(CONFIGFILE)



rule run_CBS:
    input:
        vcf = cnv_vcf_input_function('CBS')
    output:
        vcf = os.path.join(DATAPATH, "{sample_id}", "{sample_id}.CNV_calls.CBS.vcf.gz")
    threads: get_tool_resource('CBS', 'threads')
    resources:
        runtime=get_tool_resource('CBS', 'runtime'),
        mem_mb=get_tool_resource('CBS', 'memory'),
        partition=get_tool_resource('CBS', 'partition')
    params:
        # SDundo = config['settings']['CBS']['SDundo'],
        # filter=get_tool_filter_settings('CBS'),
        settings=config['settings']['CBS']
    log:
        err=os.path.join(LOGPATH, "CBS", "{sample_id}", "error.log"),
        out=os.path.join(LOGPATH, "CBS", "{sample_id}", "out.log")
    conda:
        importlib.resources.files(STEM_CNV_CHECK).joinpath("envs", "general-R.yaml")
    script:
        '../scripts/run_CBS_DNAcopy.R'
    # shell:
    #     "Rscript {SNAKEDIR}/scripts/run_CBS_DNAcopy.R -s {params.SDundo} {input.tsv} {output.tsv} {CONFIGFILE} {SAMPLETABLE} > {log.out} 2> {log.err}"




def get_processed_ref_data(wildcards):
    sample_id, ref_id = get_ref_id(wildcards)
    return os.path.join(DATAPATH, f"{ref_id}", f"{ref_id}.combined-cnv-calls.tsv") if ref_id else []


rule run_process_CNV_calls:
    input:
        cnv_calls = expand(
            os.path.join(DATAPATH, "{{sample_id}}", "{{sample_id}}.CNV_calls.{caller}.vcf.gz"),
            caller = config['settings']['CNV.calling.tools']
        ),
        ref_data = get_processed_ref_data,
        snp_vcf = cnv_vcf_input_function('settings:CNV_processing:call_processing')
    output:
        tsv = os.path.join(DATAPATH, "{sample_id}", "{sample_id}.combined-cnv-calls.tsv"),
        vcf = os.path.join(DATAPATH, "{sample_id}", "{sample_id}.combined-cnv-calls.vcf.gz")
    threads: get_tool_resource('CNV.process', 'threads')
    resources:
        runtime=get_tool_resource('CNV.process', 'runtime'),
        mem_mb=get_tool_resource('CNV.process', 'memory'),
        partition=get_tool_resource('CNV.process', 'partition')
    params:
        settings=config['settings']['CNV_processing']
    log:
        err=os.path.join(LOGPATH, "CNV_process", "{sample_id}", "error.log"),
        # out=os.path.join(LOGPATH, "CNV_process", "{sample_id}", "out.log")
    conda:
        importlib.resources.files(STEM_CNV_CHECK).joinpath("envs","general-R.yaml")
    script:
        '../scripts/process_CNV_calls.R'
    # shell:
    #     "Rscript {SNAKEDIR}/scripts/process_CNV_calls.R {params.penncnv} {params.cbs} {DATAPATH} {wildcards.sample_id} {CONFIGFILE} {SAMPLETABLE} > {log.out} 2> {log.err}"
