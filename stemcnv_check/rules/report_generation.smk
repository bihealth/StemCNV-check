import os
from stemcnv_check.helpers import config_extract, collect_SNP_cluster_ids

def get_penncnv_log_input(wildcards):
    sample_id, ref_id, sex, ref_sex = get_ref_id(wildcards, True)
    return expand(
        os.path.join(DATAPATH, sample_id, "extra_files", "PennCNV.{chrs}.error.log"),
        chrs=["auto", "chrx"] + (["chry"] if sex == "m" else []),
    )

def get_extra_snp_input_files(wildcards):
    extra_sample_def = config['evaluation_settings']['SNP_clustering']['extra_samples']
    sample_id, ref_id = get_ref_id(wildcards)

    ids = set(collect_SNP_cluster_ids(sample_id, extra_sample_def, sample_data_full))
    if ref_id:
        ids.add(ref_id)
        
    return expand(
        [
            os.path.join(
                DATAPATH,"{ids}","{ids}.annotated-SNP-data.{filter}-filter.vcf.gz"
            )
        ],
        ids=ids,
        filter=get_tool_filter_settings("evaluation_settings:SNP_clustering:filter-settings")
    )

def get_ref_excel_input(wildcards):
    sample_id, ref_id = get_ref_id(wildcards)
    if ref_id:
        return os.path.join(DATAPATH, ref_id, f"{ref_id}.summary-stats.xlsx")
    else:
        return []

rule make_summary_table:
    input:
        gencall_stats = os.path.join(DATAPATH, "{sample_id}", "extra_files", "{sample_id}.gencall-stats.txt"),
        snp_vcf = cnv_vcf_input_function("settings:CNV_processing:call_processing"),
        cnv_vcf = os.path.join(DATAPATH, "{sample_id}", "{sample_id}.combined-cnv-calls.vcf.gz"),
        penncnv_logs = get_penncnv_log_input,
        #TODO: possibly move this to a separate rule?
        extra_snp_files = get_extra_snp_input_files,
        summary_excel_ref = get_ref_excel_input,
    output:
        xlsx = os.path.join(DATAPATH, "{sample_id}", "{sample_id}.summary-stats.xlsx"),
    threads:
        lambda wildcards: len(get_extra_snp_input_files(wildcards))
    resources:
        runtime=get_tool_resource("summary_stats", "runtime"),
        mem_mb=get_tool_resource("summary_stats", "memory"),
        partition=get_tool_resource("summary_stats", "partition"),
    params:
        config = config['evaluation_settings'],
    log:
        err=os.path.join(LOGPATH, "report", "{sample_id}", "summary-stats.error.log"),
    conda:
        "../envs/general-R.yaml"
    script:
        "../scripts/summarise_stats.R"


def get_report_sample_input(wildcards):
    sample_id, ref_id, sex, ref_sex = get_ref_id(wildcards, True)
    report_settings = config["reports"][wildcards.report]

    sample_files = expand(
        [
            os.path.join(DATAPATH, "{ids}", "{ids}.combined-cnv-calls.tsv"),
            os.path.join(DATAPATH, "{ids}", "{ids}.stats.txt"),
            os.path.join(LOGPATH, "PennCNV", "{ids}", "{chrs}.error.log"),
            os.path.join(
                DATAPATH, "{ids}", "{ids}.annotated-SNP-data.{filter}-filter.vcf.gz"
            ),
        ],
        ids=(sample_id, ref_id) if ref_id else (sample_id,),
        chrs=["auto", "chrx"] + (["chry"] if sex == "m" else []),
        filter=get_tool_filter_settings(
            f"report:{wildcards.report}:call.data.and.plots"
        ),
    )

    incl_sections = config_extract(
        ("include_sections",), report_settings, config["reports"]["__default__"]
    )
    excl_sections = config_extract(
        ("exclude_sections",), report_settings, config["reports"]["__default__"]
    )
    do_snp_clustering = (
        incl_sections == "__all__" or "SNP.dendrogram" in incl_sections
    ) and not "SNP.dendrogram" in excl_sections
    if do_snp_clustering:
        extra_sample_def = config_extract(
            ("SNP_comparison", "extra_samples"),
            report_settings,
            config["reports"]["__default__"],
        )
        ids = collect_SNP_cluster_ids(sample_id, extra_sample_def, sample_data_full)
        # VCF has both filtered & unfiltered SNPs now
        # if not config_extract(('SNP_comparison', 'ignore_filter'), report_settings, config['reports']['__default__']):
        sample_files += expand(
            [
                os.path.join(
                    DATAPATH, "{ids}", "{ids}.annotated-SNP-data.{filter}-filter.vcf.gz"
                )
            ],
            ids=ids,
            filter=get_tool_filter_settings(
                f"report:{wildcards.report}:SNP_comparison"
            ),
        )

    if wildcards.ext == "pdf":
        sample_files += [os.path.join(LOGPATH, "report", "_latex_installation_check")]

    return sample_files


rule check_latex_installation:
    output:
        os.path.join(LOGPATH, "report", "_latex_installation_check"),
    resources:
        runtime=get_tool_resource("default", "runtime"),
        mem_mb=get_tool_resource("default", "memory"),
        partition=get_tool_resource("default", "partition"),
    conda:
        "../envs/general-R.yaml"
    shell:
        """
Rscript - << 'EOF'
suppressMessages(library(tinytex))

latex_path <- tinytex_root()

if (latex_path == "") {{
    message("No LaTeX installation found. Installing TinyTex.")
    install_tinytex()
}}

system("touch {output}")

EOF
"""


rule knit_report:
    input:
        get_report_sample_input,
    output:
        report=os.path.join(DATAPATH, "{sample_id}", "{sample_id}.{report}.{ext}"),
        plots=directory(os.path.join(DATAPATH, "{sample_id}", "{report}-{ext}_images")),
        # -> this extra output is optional now!
        #summary=os.path.join(DATAPATH, "{sample_id}", "{sample_id}.summary-check.tsv")
    wildcard_constraints:
        ext="(pdf|html)",
    threads: get_tool_resource("knitr", "threads")
    # params:
    #     config = get_report_config
    resources:
        runtime=get_tool_resource("knitr", "runtime"),
        mem_mb=get_tool_resource("knitr", "memory"),
        partition=get_tool_resource("knitr", "partition"),
    log:
        err=os.path.join(LOGPATH, "report", "{sample_id}", "{report}-{ext}.error.log"),
        out=os.path.join(LOGPATH, "report", "{sample_id}", "{report}-{ext}.out.log"),
    conda:
        "../envs/general-R.yaml"
    shell:
        "Rscript {SNAKEDIR}/scripts/knit_report.R {wildcards.sample_id} {wildcards.report} {CONFIGFILE} > {log.out} 2> {log.err}"
