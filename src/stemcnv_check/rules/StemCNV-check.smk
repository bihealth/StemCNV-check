# -*- coding: utf-8 -*-

import importlib.resources
import os
from pathlib import Path
from loguru import logger as logging
import tempfile
from stemcnv_check import STEM_CNV_CHECK
from stemcnv_check.helpers import read_sample_table, get_global_file
from stemcnv_check.exceptions import SampleConstraintError, ConfigValueError

SNAKEDIR = str(importlib.resources.files(STEM_CNV_CHECK))
config["snakedir"] = SNAKEDIR

SAMPLETABLE = (
    config["sample_table"] if "sample_table" in config else "sample_table.tsv"
)  # Defined by wrapper
BASEPATH = (
    config["basedir"] if "basedir" in config else os.getcwd()
)  # Defined by wrapper
DATAPATH = config[
    "data_path"
]  # if os.path.isabs(config['data_path']) else os.path.join(BASEPATH, config['data_path'])
LOGPATH = config[
    "log_path"
]  # if os.path.isabs(config['log_path']) else os.path.join(BASEPATH, config['log_path'])
TARGET = config["target"] if "target" in config else "report"  # Defined by wrapper
IDAT_INPUT = config["raw_data_folder"]


wildcard_constraints:
    sample_id=config["wildcard_constraints"]["sample_id"],
    #sentrix_name=config['wildcard_constraints']['sentrix_name'],
    #sentrix_pos=config['wildcard_constraints']['sentrix_pos'],


# Never submit these to cluster
localrules:
    all,

sample_data_df = read_sample_table(SAMPLETABLE, str(config['column_remove_regex']))

include: "common.smk"
include: "illumina_raw_processing.smk"
include: "SNP_processing.smk"
include: "penncnv.smk"
include: "report_generation.smk"

# The Call Label section of the config can have user defined labels on top of the default ones
# The default ones should always come last, but this does not happen automatically
config["evaluation_settings"]["CNV_call_labels"] = fix_call_label_order(config["evaluation_settings"]["CNV_call_labels"])


# Rule function,s move to common? ========================================================================


def get_target_files(target=TARGET):
    # Target options: ('report', 'combined-cnv-calls', 'PennCNV', 'CBS', 'SNP-data'),
    all_samples = sample_data_df['Sample_ID']

    # complete
    if target == "complete":
        out = (
            get_target_files("report") + 
            get_target_files('summary-tables') + 
            get_target_files("combined-cnv-calls") + 
            get_target_files("collate-summary")
        )
    # Report
    elif target == "report":
        out = expand(
            os.path.join(DATAPATH, "{sample_id}", "{sample_id}.{report_filetype}"),
            sample_id=all_samples,
            report_filetype=[
                rep + "." + config["reports"][rep]["file_type"]
                for rep in config["reports"].keys()
                if rep != "_default_"
            ],
        )
    # Collated summary table
    elif target == "collate-summary":
        out = [os.path.join(
            DATAPATH,
            ((config['collate_date']+'_') if 'collate_date' in config else '') + "summary-overview." + 
            config["evaluation_settings"]["collate_output"]["file_format"]
        )]
    # Stat summary tables
    elif target == "summary-tables":
        out = expand(
            os.path.join(DATAPATH, "{sample_id}", "{sample_id}.summary-stats.xlsx"),
            sample_id=all_samples,
        )
    # Collated CNV call table
    elif target == "collate-cnv-calls":
        out = [os.path.join(
            DATAPATH,
            ((config['collate_date']+'_') if 'collate_date' in config else '') + "combined-cnv-calls." +
            config["evaluation_settings"]["collate_output"]["file_format"]
        )]
    # Target Processed-calls
    elif target == "combined-cnv-calls":
        out = expand(
            os.path.join(DATAPATH, "{sample_id}", "{sample_id}.CNV_calls.combined-annotated.vcf.gz"),
            sample_id=all_samples,
        )
    # Target PennCNV
    elif target == "PennCNV":
        out = expand(
            os.path.join(DATAPATH, "{sample_id}", "{sample_id}.CNV_calls.PennCNV.vcf.gz"),
            sample_id=all_samples,
        )
    # Target CBS
    elif target == "CBS":
        out = expand(
            os.path.join(DATAPATH, "{sample_id}", "{sample_id}.CNV_calls.CBS.vcf.gz"),
            sample_id=all_samples,
        )
    # Target SNP-data
    elif target == "SNP-data":
        out = expand(
            os.path.join(
                DATAPATH,
                "{sample_id}",
                "{sample_id}.annotated-SNP-data.{filter}-filter.vcf.gz",
            ),
            sample_id=all_samples,
            filter=config["settings"]["default_probe_filter_set"],
        )
    # Target gtc-data
    elif target == "gtc-data":
        sentrix_data = sample_data_df[['Array_Name', 'Chip_Name']].drop_duplicates()
        sentrix_data['genome'] = sentrix_data.apply(
            lambda x: config['array_definition'][x['Array_Name']]['genome_version'],
            axis=1
        )
        out = [
            os.path.join(DATAPATH, f"gtc_{genome}", f"{sentrix_name}", "_done") for 
            genome, sentrix_name in zip(sentrix_data['genome'], sentrix_data['Chip_Name'])
        ]       
    else:
        raise ValueError('Invalid target value: "{}"'.format(target))

    return out


# Rules ========================================================================


rule all:
    input:
        get_target_files(),


rule run_SNV_analysis:
    input:
        snp_vcf = snp_vcf_input_function("SNV_analysis"), 
        ref_snp_vcf = get_ref_input_function(f'annotated-SNP-data.{get_tool_filter_settings("SNV_analysis")}-filter.vcf.gz'),
        extra_snp_files = get_extra_snp_input_files,
    output:
        xlsx=os.path.join(DATAPATH,"{sample_id}","{sample_id}.SNV-analysis.xlsx"),
    threads:
        1 if config['is_wsl'] else get_tool_resource("SNV_analysis","threads"),
    resources:
        runtime=get_tool_resource("SNV_analysis","runtime"),
        mem_mb=get_tool_resource("SNV_analysis","memory"),
        partition=get_tool_resource("SNV_analysis","partition"),
    params:
        config=config['settings']['SNV_analysis'],
        gtf_file= lambda wildcards: get_global_file(
            'gtf',get_static_input('genome_version')(wildcards),config['global_settings'],config['cache_path']
            ),
        ginfo_file= lambda wildcards: get_global_file(
            'genome_info',get_static_input('genome_version')(wildcards),config['global_settings'],config['cache_path']
        ),
    log:
        err=os.path.join(LOGPATH,"SNV_analysis","{sample_id}","SNV_analysis.error.log"),
    conda:
        "../envs/general-R.yaml"
    script:
        "../scripts/SNV_analysis.R"


rule run_CBS:
    input:
        vcf=snp_vcf_input_function("CBS"),
    output:
        vcf=os.path.join(DATAPATH, "{sample_id}", "{sample_id}.CNV_calls.CBS.vcf.gz"),
    threads: get_tool_resource("CBS", "threads")
    resources:
        runtime=get_tool_resource("CBS", "runtime"),
        mem_mb=get_tool_resource("CBS", "memory"),
        partition=get_tool_resource("CBS", "partition"),
    params:
        # Ensure rerun on changes to settings or sample meta data
        settings=config["settings"]["CBS"],
        sample_sex=lambda wildcards: get_sample_info(wildcards)['Sex'],
    log:
        err=os.path.join(LOGPATH, "CBS", "{sample_id}", "error.log"),
        out=os.path.join(LOGPATH, "CBS", "{sample_id}", "out.log"),
    conda:
        "../envs/general-R.yaml"
    script:
        "../scripts/run_CBS_DNAcopy.R"


rule run_process_CNV_calls:
    input:
        cnv_calls=expand(
            os.path.join(
                DATAPATH, "{{sample_id}}", "{{sample_id}}.CNV_calls.{caller}.vcf.gz"
            ),
            caller=config["settings"]["CNV.calling.tools"],
        ),
        ref_data=get_ref_input_function('CNV_calls.combined-annotated.vcf.gz'),
        snp_vcf=snp_vcf_input_function("settings:CNV_processing:call_processing"),
    output:
        vcf=os.path.join(
            DATAPATH, "{sample_id}", "{sample_id}.CNV_calls.combined-annotated.vcf.gz"
        ),
    threads: get_tool_resource("CNV.process", "threads")
    resources:
        runtime=get_tool_resource("CNV.process", "runtime"),
        mem_mb=get_tool_resource("CNV.process", "memory"),
        partition=get_tool_resource("CNV.process", "partition"),
    params:
        eval_settings=config["evaluation_settings"]["CNV_call_labels"],
        settings=config["settings"]["CNV_processing"],
        sample_roi=lambda wildcards: get_sample_info(wildcards)['Regions_of_Interest'],
        gtf_file=lambda wildcards: get_global_file(
            'gtf', get_static_input('genome_version')(wildcards), config['global_settings'], config['cache_path']
        ),
        ginfo_file=lambda wildcards: get_global_file(
            'genome_info', get_static_input('genome_version')(wildcards), config['global_settings'], config['cache_path']
        ),
        dosage_file=lambda wildcards: get_global_file(
            'dosage_scores',get_static_input('genome_version')(wildcards),config['global_settings'],config['cache_path']
        ),
        array_gaps_file=get_static_input('array_gaps_file'),
        array_density_file=get_static_input('array_density_file'),
    log:
        err=os.path.join(LOGPATH, "CNV_process", "{sample_id}", "error.log"),
    conda:
        "../envs/general-R.yaml"
    script:
        "../scripts/process_CNV_calls.R"


rule collate_cnv_calls:
    input:
        cnv_calls=expand(
            os.path.join(
                DATAPATH, "{sample_id}", "{sample_id}.CNV_calls.combined-annotated.vcf.gz"
            ),
            sample_id=sample_data_df['Sample_ID']
        ),
    output:
        os.path.join(
            DATAPATH,
            ((config['collate_date']+'_') if 'collate_date' in config else '') + "combined-cnv-calls." +
            config["evaluation_settings"]["collate_output"]["file_format"]
        ),
    threads: get_tool_resource("collate-cnv-calls", "threads")
    resources:
        runtime=get_tool_resource("collate-cnv-calls", "runtime"),
        mem_mb=get_tool_resource("collate-cnv-calls", "memory"),
        partition=get_tool_resource("collate-cnv-calls", "partition"),
    log:
        err=os.path.join(LOGPATH, "collate-cnv-calls", "error.log"),
    params:
        output_filters = config['evaluation_settings']['collate_output']['cnv_collate_call_selection'],
        output_format = config['evaluation_settings']['collate_output']['file_format'],
    conda:
        "../envs/general-R.yaml"
    script:
        "../scripts/collate_cnv_calls.R"


rule collate_summary_overview:
    input:
        summary_files=expand(
            os.path.join(DATAPATH, "{sample_id}", "{sample_id}.summary-stats.xlsx"),
            sample_id=sample_data_df['Sample_ID']
        ),
    output:
        os.path.join(
            DATAPATH,
            ((config['collate_date']+'_') if 'collate_date' in config else '') + "summary-overview." +
            config["evaluation_settings"]["collate_output"]["file_format"]
        ),
    threads: get_tool_resource("collate-summary", "threads")
    resources:
        runtime=get_tool_resource("collate-summary", "runtime"),
        mem_mb=get_tool_resource("collate-summary", "memory"),
        partition=get_tool_resource("collate-summary", "partition"),
    # log:
    #     err=os.path.join(LOGPATH, "summary-overview", "error.log"),
    params:
        sampletable=sample_data_df,
        extra_cols=config['evaluation_settings']['collate_output']['summary_extra_sampletable_cols'],
        output_format= config['evaluation_settings']['collate_output']['file_format'],
    run:
        import pandas as pd
        def prep_summary(file):
            df = pd.read_excel(file, sheet_name='summary_stats').T
            df = (
                df.rename(columns=df.iloc[0]).
                drop(df.index[0]).
                head(2).
                reset_index().
                melt(id_vars=('sample_id', 'index')).
                assign(index = lambda x: x['index'].str.replace('sample_', '').replace('value', '')).
                pivot(index='sample_id', columns=['variable', 'index'], values='value')
                )
            df.columns = ['{a}{sep}{b}'.format(a=a, b=b, sep="_" if b else "") for a, b in df.columns]
            return df
        
        summary = (pd.concat([prep_summary(file) for file in input.summary_files]).
                   drop(columns=['SNPs_post_filter_eval']))
        # Add extra columns from sampletable
        sample_info = params.sampletable[['Sample_ID'] + params.extra_cols].set_index('Sample_ID')
        summary = sample_info.join(summary, validate='1:1')
        
        if params.output_format == 'xlsx':
            summary.to_excel(output[0], index=True, index_label='Sample_ID')
        else: #tsv
            summary.to_csv(output[0], index=True, index_label='Sample_ID', sep='\t')