import os

localrules:
    relink_gencall,

rule run_gencall:
    input:
        bpm=get_static_input("bpm_manifest_file"),
        egt=get_static_input("egt_cluster_file"),
        idat_path=os.path.join(IDAT_INPUT, "{sentrix_name}"),
    output:
        os.path.join(DATAPATH, "gtc_{genome}", "{sentrix_name}", "_done"),
    threads: get_tool_resource("GenCall", "threads")
    resources:
        runtime=get_tool_resource("GenCall", "runtime"),
        mem_mb=get_tool_resource("GenCall", "memory"),
        partition=get_tool_resource("GenCall", "partition"),
    params:
        options='--gender-estimate-call-rate-threshold -0.1',
        outpath=lambda wildcards: fix_container_path(
            os.path.join(DATAPATH, f"gtc_{wildcards.genome}", wildcards.sentrix_name), "data"
        ),
        bpm=lambda wildcards: fix_container_path(
            get_static_input("bpm_manifest_file")(wildcards),
            # The get_sample_info function requires the sample_id wildcard, also need to deal with likely duplicated Chip_Names
            set(sample_data_df.loc[sample_data_df['Chip_Name'] == wildcards.sentrix_name]['Array_Name'].values).pop()
        ),
        egt=lambda wildcards: fix_container_path(
            get_static_input("egt_cluster_file")(wildcards),
            set(sample_data_df.loc[sample_data_df['Chip_Name'] == wildcards.sentrix_name]['Array_Name'].values).pop()
        ),
        idat_path=lambda wildcards: fix_container_path(
            os.path.join(IDAT_INPUT, wildcards.sentrix_name), "rawdata"
        ),
        logerr=lambda wildcards: fix_container_path(
            os.path.join(LOGPATH, "GenCall", wildcards.sentrix_name, f"{wildcards.genome}.error.log"),
            "logs",
        ),
        logout=lambda wildcards: fix_container_path(
            os.path.join(LOGPATH, "GenCall", wildcards.sentrix_name, f"{wildcards.genome}.out.log"), "logs"
        ),
    log:
        err=os.path.join(LOGPATH, "GenCall", "{sentrix_name}", "{genome}.error.log"),
        out=os.path.join(LOGPATH, "GenCall", "{sentrix_name}", "{genome}.out.log"),
    container:
        # Not sure if we need to use a specific version here
        "docker://us.gcr.io/broad-gotc-prod/illumina-iaap-autocall:1.0.2-1.1.0-1629910298"
    shell:
        '/usr/gitc/iaap/iaap-cli/iaap-cli gencall '
        '"{params.bpm}" "{params.egt}" "{params.outpath}" '
        '--idat-folder "{params.idat_path}" '
        '--output-gtc {params.options} '
        '-t {threads} '
        '> {params.logout} '
        '2> {params.logerr} '
        '&& [ $(ls {params.outpath}/*.gtc -l | wc -l) -ge 1 ] '
        '&& touch {params.outpath}/_done || exit 1'


# The iaap-cli will *always* generate filenames derived from the sentrix name & pos
def get_chip(wildcards, outtype="dir_path"):
    """Get the chip name from a sample_id
    Values for outtype: 'dirpath' | 'file'"""
    chip_name = get_sample_info(wildcards)['Chip_Name']
    chip_pos = get_sample_info(wildcards)['Chip_Pos']
    genome = get_static_input("genome_version")(wildcards)
    if outtype == "dir_path":
        return os.path.join(DATAPATH, f"gtc_{genome}", chip_name)
    elif outtype == "file":
        return os.path.join(chip_name, chip_name + "_" + chip_pos + ".gtc")


# Note:
# the *.gtc output files from gencall will *always* match the idat_file names
# --> maybe better to switch from Chip_Name & Chip_Pos to Folder_Name & IDAT_Name ?!
rule relink_gencall:
    input:
        lambda wildcards: os.path.join(get_chip(wildcards), "_done"),
    output:
        os.path.join(DATAPATH, "{sample_id}", "{sample_id}.gencall.{genome}.gtc"),
    params:
        gtc_link_path=lambda wildcards: os.path.join(
            "..", f"gtc_{wildcards.genome}", get_chip(wildcards, outtype="file")
        ),
    shell:
        'ln -s "{params.gtc_link_path}" "{output}"'


rule run_gtc2vcf_vcf:
    input:
        bpm=get_static_input("bpm_manifest_file"),
        egt=get_static_input("egt_cluster_file"),
        genome=get_static_input('fasta'),
        gtc=lambda wildcards: os.path.join(
            DATAPATH, f"{wildcards.sample_id}", f"{wildcards.sample_id}.gencall.{get_static_input("genome_version")(wildcards)}.gtc"
        ),
    output:
        vcf = pipe(os.path.join(DATAPATH,"{sample_id}","{sample_id}.unprocessed.vcf")),
        stats = os.path.join(DATAPATH, "{sample_id}", "extra_files", "{sample_id}.gencall-stats.txt")
    threads: get_tool_resource("gtc2vcf", "threads")
    resources:
        runtime=get_tool_resource("gtc2vcf", "runtime"),
        mem_mb=get_tool_resource("gtc2vcf", "memory"),
        partition=get_tool_resource("gtc2vcf", "partition"),
    params:
        csv=lambda wildcards: (
            '--csv "{}"'.format(get_static_input("csv_manifest_file")(wildcards))
            if get_static_input("csv_manifest_file")(wildcards)
            else ""
        ),
    log:
        vcf=os.path.join(LOGPATH, "gtc2vcf", "{sample_id}", "bcftools.gtc2vcf.log"),
        sort=os.path.join(LOGPATH, "gtc2vcf", "{sample_id}", "bcftools.sort.log"),
        norm=os.path.join(LOGPATH, "gtc2vcf", "{sample_id}", "bcftools.norm.log"),
    conda:
        "../envs/gtc2vcf.yaml"
    shell:
        'bcftools plugin gtc2vcf -O v '
        '--bpm "{input.bpm}" {params.csv} --egt "{input.egt}" '
        '--fasta-ref "{input.genome}" --extra {output.stats} '
        '{input.gtc} 2> {log.vcf} | '
        'bcftools sort 2> {log.sort} | '
        'bcftools norm -m- --multi-overlaps . -o {output.vcf} 2> {log.norm}'
