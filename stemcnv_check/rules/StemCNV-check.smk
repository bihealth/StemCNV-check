# -*- coding: utf-8 -*-

import importlib.resources
import os
from pathlib import Path
from loguru import logger as logging
import tempfile
import ruamel.yaml as ruamel_yaml
from stemcnv_check import STEM_CNV_CHECK
from stemcnv_check.helpers import read_sample_table
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
    f, CONFIGFILE = tempfile.mkstemp(suffix=".yaml", text=True)
    os.close(f)
    removetempconfig = True


config["snakedir"] = SNAKEDIR
# config["genome_version"] = config["genome_version"][-2:]

yaml = ruamel_yaml.YAML(typ="safe")
with open(CONFIGFILE, "w") as yamlout:
    yaml.dump(config, yamlout)


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


sample_data = read_sample_table(SAMPLETABLE)
sample_data_full = read_sample_table(SAMPLETABLE, with_opt=True)


include: "common.smk"
include: "illumina_raw_processing.smk"
include: "SNP_processing.smk"
include: "penncnv.smk"
include: "report_generation.smk"


# Rule function,s move to common? ========================================================================


def get_target_files(target=TARGET):
    # Target options: ('report', 'combined-cnv-calls', 'PennCNV', 'CBS', 'SNP-data'),
    all_samples = [sample_id for sample_id, _, _, _, _ in sample_data]

    # complete
    if target == "complete":
        out = get_target_files("report") + get_target_files("combined-cnv-calls")

    # Report
    elif target == "report":
        out = expand(
            os.path.join(DATAPATH, "{sample_id}", "{sample_id}.{report_filetype}"),
            sample_id=all_samples,
            report_filetype=[
                rep + "." + config["reports"][rep]["file_type"]
                for rep in config["reports"].keys()
                if rep != "__default__"
            ],
        )  # + \
        # expand(os.path.join(DATAPATH,"{sample_id}","{sample_id}.summary-check.tsv"),
        #                sample_id = all_samples)
    # Target Processed-calls
    elif target == "combined-cnv-calls":
        out = expand(
            [
                # os.path.join(
                #     DATAPATH, "{sample_id}", "{sample_id}.combined-cnv-calls.tsv"
                # ),
                os.path.join(
                    DATAPATH, "{sample_id}", "{sample_id}.combined-cnv-calls.vcf.gz"
                ),
            ],
            sample_id=all_samples,
        )
    # Target PennCNV
    elif target == "PennCNV":
        out = expand(
            os.path.join(
                DATAPATH, "{sample_id}", "{sample_id}.CNV_calls.penncnv.vcf.gz"
            ),
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
        # TODO: update this
        out = expand(
            os.path.join(
                DATAPATH,
                "{sample_id}",
                "{sample_id}.annotated-SNP-data.{filter}-filter.vcf.gz",
            ),
            # os.path.join(DATAPATH, "{sample_id}", "{sample_id}.annotated-SNP-data.{filter}-filter.vcf.gz"),
            sample_id=all_samples,
            filter=config["settings"]["default-filter-set"],
        )
        
    else:
        raise ValueError('Invalid target value: "{}"'.format(target))

    return out


# Rules ========================================================================


rule all:
    input:
        get_target_files(),
    run:
        if removetempconfig:
            os.remove(CONFIGFILE)


rule run_CBS:
    input:
        vcf=cnv_vcf_input_function("CBS"),
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
        sex_info=lambda wildcards: get_ref_id(wildcards,True),
    log:
        err=os.path.join(LOGPATH, "CBS", "{sample_id}", "error.log"),
        out=os.path.join(LOGPATH, "CBS", "{sample_id}", "out.log"),
    conda:
        "../envs/general-R.yaml"
    script:
        "../scripts/run_CBS_DNAcopy.R"


# shell:
#     "Rscript {SNAKEDIR}/scripts/run_CBS_DNAcopy.R -s {params.SDundo} {input.tsv} {output.tsv} {CONFIGFILE} {SAMPLETABLE} > {log.out} 2> {log.err}"


def get_processed_ref_data(wildcards):
    sample_id, ref_id = get_ref_id(wildcards)
    return (
        os.path.join(DATAPATH, f"{ref_id}", f"{ref_id}.combined-cnv-calls.vcf.gz")
        if ref_id
        else []
    )


rule run_process_CNV_calls:
    input:
        cnv_calls=expand(
            os.path.join(
                DATAPATH, "{{sample_id}}", "{{sample_id}}.CNV_calls.{caller}.vcf.gz"
            ),
            caller=config["settings"]["CNV.calling.tools"],
        ),
        ref_data=get_processed_ref_data,
        snp_vcf=cnv_vcf_input_function("settings:CNV_processing:call_processing"),
    output:
        # tsv=os.path.join(DATAPATH, "{sample_id}", "{sample_id}.combined-cnv-calls.tsv"),
        vcf=os.path.join(
            DATAPATH, "{sample_id}", "{sample_id}.combined-cnv-calls.vcf.gz"
        ),
    threads: get_tool_resource("CNV.process", "threads")
    resources:
        runtime=get_tool_resource("CNV.process", "runtime"),
        mem_mb=get_tool_resource("CNV.process", "memory"),
        partition=get_tool_resource("CNV.process", "partition"),
    params:
        settings=config["settings"]["CNV_processing"],
    log:
        err=os.path.join(LOGPATH, "CNV_process", "{sample_id}", "error.log"),
        # out=os.path.join(LOGPATH, "CNV_process", "{sample_id}", "out.log")
    conda:
        "../envs/general-R.yaml"
    script:
        "../scripts/process_CNV_calls.R"
