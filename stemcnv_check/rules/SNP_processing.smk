import os
import importlib.resources
from stemcnv_check import STEM_CNV_CHECK
from stemcnv_check.helpers import get_mehari_db_file


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
        sample_sex=lambda wildcards: get_ref_id(wildcards, True)[2]
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
        genomeversion=(
            "grch38" if config["genome_version"] in ("hg38", "GRCh38") else "grch37"
        ),
        mehari_db_path=get_mehari_db_file(
            config['global_settings']['mehari_transcript_db'],
            config['cache_path'],
            config["genome_version"]
        ),
    conda:
        "../envs/snp-annotation.yaml"
    shell:
        "mehari annotate seqvars "
        "--path-input-vcf {input.vcf} "
        "--path-output-vcf {output} "
        "--transcripts {params.mehari_db_path} "
        "--genome-release {params.genomeversion} "
        "--pick-transcript-mode first "
        "--report-most-severe-consequence-by allele "
        #"--hgnc <path> "
        " > {log.out} 2> {log.err}"
