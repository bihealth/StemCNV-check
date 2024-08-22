
rule run_gencall:
    input:
        bpm=config["static_data"]["bpm_manifest_file"],
        egt=config["static_data"]["egt_cluster_file"],
        idat_path=os.path.join(IDAT_INPUT, "{sentrix_name}"),
    output:
        os.path.join(DATAPATH, "gtc", "{sentrix_name}", "_done"),
    threads: get_tool_resource("GenCall", "threads")
    resources:
        runtime=get_tool_resource("GenCall", "runtime"),
        mem_mb=get_tool_resource("GenCall", "memory"),
        partition=get_tool_resource("GenCall", "partition"),
    params:
        options=get_tool_resource("GenCall", "cmd-line-params"),
        outpath=lambda wildcards: fix_container_path(
            os.path.join(DATAPATH, "gtc", wildcards.sentrix_name), "data"
        ),
        bpm=fix_container_path(config["static_data"]["bpm_manifest_file"], "static"),
        egt=fix_container_path(config["static_data"]["egt_cluster_file"], "static"),
        idat_path=lambda wildcards: fix_container_path(
            os.path.join(IDAT_INPUT, wildcards.sentrix_name), "rawdata"
        ),
        logerr=lambda wildcards: fix_container_path(
            os.path.join(LOGPATH, "GenCall", wildcards.sentrix_name, "error.log"),
            "logs",
        ),
        logout=lambda wildcards: fix_container_path(
            os.path.join(LOGPATH, "GenCall", wildcards.sentrix_name, "out.log"), "logs"
        ),
    log:
        err=os.path.join(LOGPATH, "GenCall", "{sentrix_name}", "error.log"),
        out=os.path.join(LOGPATH, "GenCall", "{sentrix_name}", "out.log"),
    container:
        # Not sure if we need to use a specific version here
        "docker://us.gcr.io/broad-gotc-prod/illumina-iaap-autocall:1.0.2-1.1.0-1629910298"
    shell:
        '/usr/gitc/iaap/iaap-cli/iaap-cli gencall "{params.bpm}" "{params.egt}" "{params.outpath}" --idat-folder "{params.idat_path}" --output-gtc {params.options} -t {threads} > {params.logout} 2> {params.logerr} && [ $(ls {params.outpath}/*.gtc -l | wc -l) -ge 1 ] && touch {params.outpath}/_done || exit 1'


# The iaap-cli will *always* generate filenames derived from the sentrix name & pos
def get_chip(wildcards, outtype="dir_path"):
    """Get the chip name from a sample_id
    Values for outtype: 'dirpath' | 'file'"""
    chip_name, chip_pos = [
        (n, p) for sid, n, p, _, _ in sample_data if sid == wildcards.sample_id
    ][0]
    if outtype == "dir_path":
        return os.path.join(DATAPATH, "gtc", chip_name)
    elif outtype == "file":
        return os.path.join(chip_name, chip_name + "_" + chip_pos + ".gtc")


# Note:
# the *.gtc output files from gencall will *always* match the idat_file names
# --> maybe better to switch from Chip_Name & Chip_Pos to Folder_Name & IDAT_Name ?!
rule relink_gencall:
    input:
        lambda wildcards: os.path.join(get_chip(wildcards), "_done"),
    output:
        os.path.join(DATAPATH, "{sample_id}", "{sample_id}.gencall.gtc"),
    params:
        gtc_link_path=lambda wildcards: os.path.join(
            "..", "gtc", get_chip(wildcards, outtype="file")
        ),
    shell:
        'ln -s "{params.gtc_link_path}" "{output}"'


# TODO: input functions to get correct fasta (& later correct manifest files)
rule run_gtc2vcf_vcf:
    input:
        bpm=config["static_data"]["bpm_manifest_file"],
        egt=config["static_data"]["egt_cluster_file"],
        genome=get_genome_fasta,
        gtc=os.path.join(DATAPATH, "{sample_id}", "{sample_id}.gencall.gtc"),
    output:
        vcf=pipe(os.path.join(DATAPATH, "{sample_id}", "{sample_id}.unprocessed.vcf")),
        metatxt=os.path.join(DATAPATH, "{sample_id}", "{sample_id}.stats.txt"),
    threads: get_tool_resource("gtc2vcf", "threads")
    resources:
        runtime=get_tool_resource("gtc2vcf", "runtime"),
        mem_mb=get_tool_resource("gtc2vcf", "memory"),
        partition=get_tool_resource("gtc2vcf", "partition"),
    params:
        options=get_tool_resource("gtc2vcf", "cmd-line-params"),
        csv=(
            '--csv "{}"'.format(config["static_data"]["csv_manifest_file"])
            if config["static_data"]["csv_manifest_file"]
            else ""
        ),
    log:
        err=os.path.join(LOGPATH, "gtc2vcf", "{sample_id}", "vcf.error.log"),
        #out=os.path.join(LOGPATH, "gtc2vcf", "{sample_id}", "vcf.out.log"),
    conda:
        "../envs/gtc2vcf.yaml"
    shell:
        'bcftools plugin gtc2vcf {params.options} -O v --bpm "{input.bpm}" {params.csv} --egt "{input.egt}" --fasta-ref "{input.genome}" --extra {output.metatxt} {input.gtc} 2> {log.err} | bcftools sort -o {output.vcf} 2>> {log.err}'
