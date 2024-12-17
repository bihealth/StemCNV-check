import os
import importlib.resources
from stemcnv_check import STEM_CNV_CHECK
from stemcnv_check.helpers import get_global_file


rule filter_snp_vcf:
    input:
        os.path.join(DATAPATH, "{sample_id}", "{sample_id}.unprocessed.vcf"),
    output:
        pipe(os.path.join(
            DATAPATH,
            "{sample_id}",
            "{sample_id}.processed-SNP-data.{filter}-filter.vcf",
        )),
    threads: get_tool_resource("filter_snp_vcf", "threads")
    resources:
        runtime=get_tool_resource("filter_snp_vcf", "runtime"),
        mem_mb=get_tool_resource("filter_snp_vcf", "memory"),
        partition=get_tool_resource("filter_snp_vcf", "partition"),
    log:
        err=os.path.join(LOGPATH, "filter_snp_vcf", "{sample_id}", "{filter}.error.log"),
        #out=os.path.join(LOGPATH, "filter_snp_vcf", "{sample_id}", "{filter}.out.log")
    params:
        filter=lambda wildcards: config['settings']['probe-filter-sets'][wildcards.filter],
        sample_sex=lambda wildcards: get_sample_info(wildcards)['Sex'],
        genome_version=lambda wildcards: get_static_input('genome_version')(wildcards)
    conda:
        "../envs/python-vcf.yaml"
    script:
        "../scripts/filter_snp_vcf.py"

rule mehari_annotate_snp_vcf:
    input:
        vcf=os.path.join(
            DATAPATH,
            "{sample_id}",
            "{sample_id}.processed-SNP-data.{filter}-filter.vcf",
        ),
    output:
        os.path.join(
            DATAPATH,
            "{sample_id}",
            "{sample_id}.annotated-SNP-data.{filter}-filter.vcf.gz",
        ),
    resources:
        threads=get_tool_resource("mehari", "threads"),
        runtime=get_tool_resource("mehari", "runtime"),
        mem_mb=get_tool_resource("mehari", "memory"),
        partition=get_tool_resource("mehari", "partition"),
    log:
        err=os.path.join(
            LOGPATH, "mehari_annotate", "{sample_id}", "snp_vcf.{filter}.error.log"
        ),
        out=os.path.join(LOGPATH, "mehari_annotate", "{sample_id}", "snp_vcf.{filter}.out.log"),
    params:
        genomeversion=lambda wildcards: (
            "grch38" if get_static_input('genome_version')(wildcards) in ("hg38", "GRCh38") else "grch37"
        ),
        mehari_db_path=lambda wildcards: get_static_input('mehari_txdb')(wildcards)
    conda:
        "../envs/snp-annotation.yaml"
    shell:
        "mehari annotate seqvars "
        "--path-input-vcf {input.vcf} "
        "--path-output-vcf {output} "
        "--transcripts {params.mehari_db_path} "
        "--genome-release {params.genomeversion} "
        # need to select transcript pick order, for pick-mode to work
        "--pick-transcript mane-select  "
        "--pick-transcript mane-plus-clinical "
        "--pick-transcript length "
        "--pick-transcript-mode first "
        "--report-most-severe-consequence-by allele "
        #"--hgnc <path> "
        " > {log.out} 2> {log.err}"
