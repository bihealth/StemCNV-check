import os
import importlib.resources
from stemcnv_check import STEM_CNV_CHECK


rule filter_snp_vcf:
    input:
        os.path.join(DATAPATH, "{sample_id}", "{sample_id}.unprocessed.vcf"),
    output:
        os.path.join(
            DATAPATH,
            "{sample_id}",
            "{sample_id}.processed-SNP-data.{filter}-filter.vcf",
        ),
    threads: get_tool_resource("filter_snp_vcf", "threads")
    resources:
        runtime=get_tool_resource("filter_snp_vcf", "runtime"),
        mem_mb=get_tool_resource("filter_snp_vcf", "memory"),
        partition=get_tool_resource("filter_snp_vcf", "partition"),
    log:
        err=os.path.join(LOGPATH, "filter_snp_vcf", "{sample_id}", "{filter}.error.log"),
        #out=os.path.join(LOGPATH, "filter_snp_vcf", "{sample_id}", "{filter}.out.log")
    # conda:
    #     importlib.resources.files(STEM_CNV_CHECK).joinpath("envs","general-R.yaml")
    # script:
    #     '../scripts/filter_snp_vcf.R'
    conda:
        "../envs/python-vcf.yaml"
    script:
        "../scripts/filter_snp_vcf.py"


rule annotate_snp_vcf:
    input:
        vcf=os.path.join(
            DATAPATH,
            "{sample_id}",
            "{sample_id}.processed-SNP-data.{filter}-filter.vcf",
        ),
        genomefasta=get_genome_fasta,
    output:
        os.path.join(
            DATAPATH,
            "{sample_id}",
            "{sample_id}.annotated-SNP-data.{filter}-filter.vcf.gz",
        ),
    threads: get_tool_resource("VEP", "threads")
    resources:
        threads=get_tool_resource("VEP", "threads"),
        runtime=get_tool_resource("VEP", "runtime"),
        mem_mb=get_tool_resource("VEP", "memory"),
        partition=get_tool_resource("VEP", "partition"),
    log:
        err=os.path.join(
            LOGPATH, "annotate_snp_vcf", "{sample_id}", "{filter}.error.log"
        ),
        out=os.path.join(LOGPATH, "annotate_snp_vcf", "{sample_id}", "{filter}.out.log"),
    params:
        genomeversion=(
            "GRCh38" if config["genome_version"] in ("hg38", "GRCh38") else "GRCh37"
        ),
        vep_cache_path=config["use_vep_cache"],
    conda:
        importlib.resources.files(STEM_CNV_CHECK).joinpath(
            "envs", "vep-annotation.yaml"
        )
    shell:
        "vep --verbose "
        "--fasta {input.genomefasta} "

        "--input_file {input.vcf} "
        "--output_file {output} "
        "--compress_output bgzip "
        "--format vcf "
        "--vcf "
        "--force_overwrite "
        "--no_stats "
        "--warning_file {log.err} "
        "--skipped_variants_file {log.out} "
        "--assembly {params.genomeversion} "
        "--fork {resources.threads} "
        "--cache "
        "--dir_cache {params.vep_cache_path} "

        "--total_length "

        "--gencode_basic "
        "--symbol "
        "--terms SO "
        "--hgvs "
        "--pick "


        "--check_existing --no_check_alleles "
        "--af "
        "--pubmed "

        '--fields "Gene,SYMBOL,STRAND,Consequence,cDNA_position,CDS_position,Protein_position,HGVSc,HGVSp,Existing_variation,CLIN_SIG,SOMATIC,PHENO,AF,PUBMED" '

        " >> {log.out} 2>> {log.err}"
        # only needed for HGVS
        # TODO: maye omitting this helps w/ pipe issues?
        #might help with pipe auto-detection issues
        ## Annotation options
        # Gene & protein annotation
        #> limit to gencode transcripts
        #> add gene symbol
        # how to write/format/annotate the consequence
        #> add HGVS nomenclature (protein changes)
        #> pick the most severe consequence (& gene) per variant
        # will query databases for existing annotation at the same position
        # includes existing annotations from ClinVar, COSMIC etc
        #  global allele frequency (AF) from 1000 Genomes Phase 3 data 
        # Select content of CSQ:
