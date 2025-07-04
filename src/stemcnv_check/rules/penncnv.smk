import os

# Never submit these to cluster
localrules:
    prep_PennCNV_sexfile,


rule prep_PennCNV_sexfile:
    input:
        snp_vcf_input_function("PennCNV"),
    output:
        temp(os.path.join(DATAPATH, "{sample_id}", "{sample_id}.penncnv.sexfile.txt")),
    params:
        sample_sex=lambda wildcards: get_sample_info(wildcards)['Sex'],
        sample_docker_path=lambda wildcards: fix_container_path(
            os.path.join(
                DATAPATH,
                wildcards.sample_id,
                f"{wildcards.sample_id}.penncnv.input.tsv",
            ),
            "data",
        ),
    shell:
        'echo -e "{params.sample_docker_path}\t{params.sample_sex}" > {output}; '


# - extract tsv SNP file from vcf
rule prep_PennCNV_input:
    input:
        vcf=snp_vcf_input_function("PennCNV"),
    output:
        tsv=temp(os.path.join(DATAPATH, "{sample_id}", "{sample_id}.penncnv.input.tsv")),
    log:
        os.path.join(LOGPATH, "PennCNV", "{sample_id}", "input.log"),
    resources:
        runtime=get_tool_resource("PennCNV", "runtime"),
        mem_mb=get_tool_resource("PennCNV", "memory"),
        partition=get_tool_resource("PennCNV", "partition"),
    conda:
        "../envs/vembrane.yaml"
    params:
        filter=get_tool_filter_settings("PennCNV"),
    shell:
        "vembrane filter '\"PASS\" in FILTER' {input.vcf} 2> {log} | "
        "vembrane table --header 'Name, Chr, Position, B Allele Freq, Log R Ratio'"
        ' \'ID, CHROM, POS, FORMAT["BAF"][SAMPLE], FORMAT["LRR"][SAMPLE]\' > {output.tsv} 2>> {log}'
# # Alternative:
#     wrapper:
#         "v3.14.1/bio/vembrane/table"


# - run auto, x & Y calling
rule run_PennCNV:
    input:
        tsv=os.path.join(DATAPATH, "{sample_id}", "{sample_id}.penncnv.input.tsv"),
        sexfile=os.path.join(DATAPATH, "{sample_id}", "{sample_id}.penncnv.sexfile.txt"),
        pfb=get_static_input("penncnv_pfb_file"),
        gcmodel=get_static_input("penncnv_GCmodel_file"),
    output:
        tsv=temp(os.path.join(DATAPATH, "{sample_id}", "{sample_id}.penncnv-{chr}.tsv")),
        extra=os.path.join(DATAPATH, "{sample_id}", "extra_files", "PennCNV.{chr}.error.log"),
    threads: get_tool_resource("PennCNV", "threads")
    resources:
        runtime=get_tool_resource("PennCNV", "runtime"),
        mem_mb=get_tool_resource("PennCNV", "memory"),
        partition=get_tool_resource("PennCNV", "partition"),
    wildcard_constraints:
        chr="chrx|chry|auto",
    params:
        filter=get_tool_filter_settings("PennCNV"),
        chrom=lambda wildcards: "" if wildcards.chr == "auto" else "-" + wildcards.chr,
        #LOH calling is optional, also disabled for male sex chromosomes
        do_loh=lambda wildcards: (
            ""
            if (
                (wildcards.chr != "auto"
                and get_sample_info(wildcards)['Sex'] == "m")
                or not config["settings"]["PennCNV"]["enable_LOH_calls"]
            )
            else "-loh"
        ),
        snakedir=fix_container_path(SNAKEDIR, "snakedir"),
        tsvout=lambda wildcards: fix_container_path(
            os.path.join(
                DATAPATH,
                wildcards.sample_id,
                f"{wildcards.sample_id}.penncnv-{wildcards.chr}.tsv",
            ),
            "data",
        ),
        extraout=lambda wildcards: fix_container_path(
            os.path.join(
                DATAPATH,
                wildcards.sample_id,
                "extra_files",
                f"PennCNV.{wildcards.chr}.error.log",
            ),
            "data",
        ),
        tsvin=lambda wildcards: fix_container_path(
            os.path.join(
                DATAPATH,
                wildcards.sample_id,
                f"{wildcards.sample_id}.penncnv.input.tsv",
            ),
            "data",
        ),
        pfb=lambda wildcards: fix_container_path(
            get_static_input("penncnv_pfb_file")(wildcards), 
            get_sample_info(wildcards)["Array_Name"]
        ),
        sexfile=lambda wildcards: (
            ""
            if wildcards.chr == "auto"
            else "--sexfile "
            + str(
                fix_container_path(
                    os.path.join(
                        DATAPATH,
                        wildcards.sample_id,
                        f"{wildcards.sample_id}.penncnv.sexfile.txt",
                    ),
                    "data",
                )
            )
        ),
        gcmodel=lambda wildcards: fix_container_path(
            get_static_input("penncnv_GCmodel_file")(wildcards), 
            get_sample_info(wildcards)["Array_Name"]
        ),
        logerr=lambda wildcards: fix_container_path(
            os.path.join(
                LOGPATH, "PennCNV", wildcards.sample_id, f"{wildcards.chr}.error.log"
            ),
            "logs",
        ),
        logout=lambda wildcards: fix_container_path(
            os.path.join(
                LOGPATH, "PennCNV", wildcards.sample_id, f"{wildcards.chr}.out.log"
            ),
            "logs",
        ),
    log:
        err=os.path.join(LOGPATH, "PennCNV", "{sample_id}", "{chr}.error.log"),
        out=os.path.join(LOGPATH, "PennCNV", "{sample_id}", "{chr}.out.log"),
    container:
        "docker://genomicslab/penncnv"
    shell:
        "/home/user/PennCNV/detect_cnv.pl -test {params.do_loh} {params.chrom} -confidence "
        "-hmm {params.snakedir}/supplemental-files/hhall_loh.hmm -pfb {params.pfb} -gcmodel {params.gcmodel} "
        "{params.sexfile} {params.tsvin} -out {params.tsvout} > {params.logout} 2> {params.logerr} && "
        "cp {params.logerr} {params.extraout}"


def get_penncnv_output(wildcards, files="tsv"):
    sex = get_sample_info(wildcards)['Sex']
    chrs = ["auto", "chrx"]
    if sex == "m":
        chrs.append("chry")

    if files == "tsv":
        return expand(
            os.path.join(DATAPATH, "{sample_id}", "{sample_id}.penncnv-{chr}.tsv"),
            chr=chrs,
            sample_id=wildcards.sample_id,
        )
    elif files == "log":
        return expand(
            os.path.join(LOGPATH, "PennCNV", "{sample_id}", "{chr}.error.log"),
            chr=chrs,
            sample_id=wildcards.sample_id,
        )
    else:
        raise ValueError("Invalid file type requested: {}".format(files))


# - combine calls and write out as vcf
rule combined_PennCNV_output:
    input:
        vcf=snp_vcf_input_function("PennCNV"),
        tsvs=get_penncnv_output,
        #logs=get_penncnv_output(wildcards, 'log')
    output:
        vcf=os.path.join(
            DATAPATH, "{sample_id}", "{sample_id}.CNV_calls.PennCNV.vcf.gz"
        ),
        # FIXME (future): ideally handle this in separate
        # stats=os.path.join(DATAPATH,"{sample_id}","{sample_id}.CNV_calls.penncnv.stats.tsv")
    log:
        err=os.path.join(LOGPATH, "PennCNV", "{sample_id}", "combine.error.log"),
    resources:
        runtime=get_tool_resource("PennCNV", "runtime"),
        mem_mb=get_tool_resource("PennCNV", "memory"),
        partition=get_tool_resource("PennCNV", "partition"),
    conda:
        "../envs/general-R.yaml"
    script:
        "../scripts/combine_penncnv.R"
