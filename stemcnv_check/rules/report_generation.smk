import importlib.resources
import os
from deepdiff import DeepDiff
from pydantic.v1.utils import deep_update
from stemcnv_check.helpers import load_config, collect_SNP_cluster_ids
from stemcnv_check import __version__, STEM_CNV_CHECK

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

rule make_summary_table:
    input:
        gencall_stats = os.path.join(DATAPATH, "{sample_id}", "extra_files", "{sample_id}.gencall-stats.txt"),
        snp_vcf = cnv_vcf_input_function("settings:CNV_processing:call_processing"),
        cnv_vcf = os.path.join(DATAPATH, "{sample_id}", "{sample_id}.combined-cnv-calls.vcf.gz"),
        penncnv_logs = get_penncnv_log_input,
        #TODO: possibly move this to a separate rule?
        extra_snp_files = get_extra_snp_input_files,
        summary_excel_ref = get_ref_input_function('summary-stats.xlsx'),
    output:
        xlsx = os.path.join(DATAPATH, "{sample_id}", "{sample_id}.summary-stats.xlsx"),
    threads:
        get_tool_resource("summary_stats", "threads"),
    resources:
        runtime=get_tool_resource("summary_stats", "runtime"),
        mem_mb=get_tool_resource("summary_stats", "memory"),
        partition=get_tool_resource("summary_stats", "partition"),
    params:
        config = config['evaluation_settings']
    log:
        err=os.path.join(LOGPATH, "report", "{sample_id}", "summary-stats.error.log"),
    conda:
        "../envs/general-R.yaml"
    script:
        "../scripts/summarise_stats.R"


# def get_report_sample_input(wildcards):
#     sample_id, ref_id, sex, ref_sex = get_ref_id(wildcards, True)
#     report_settings = config["reports"][wildcards.report]
# 
#     sample_files = expand(
#         [
#             os.path.join(DATAPATH, "{ids}", "{ids}.combined-cnv-calls.tsv"),
#             os.path.join(DATAPATH, "{ids}", "{ids}.stats.txt"),
#             os.path.join(LOGPATH, "PennCNV", "{ids}", "{chrs}.error.log"),
#             os.path.join(
#                 DATAPATH, "{ids}", "{ids}.annotated-SNP-data.{filter}-filter.vcf.gz"
#             ),
#         ],
#         ids=(sample_id, ref_id) if ref_id else (sample_id,),
#         chrs=["auto", "chrx"] + (["chry"] if sex == "m" else []),
#         filter=get_tool_filter_settings(
#             f"report:{wildcards.report}:call.data.and.plots"
#         ),
#     )
# 
#     incl_sections = config_extract(
#         ("include_sections",), report_settings, config["reports"]["_default_"]
#     )
#     excl_sections = config_extract(
#         ("exclude_sections",), report_settings, config["reports"]["_default_"]
#     )
#     do_snp_clustering = (
#         incl_sections == "__all__" or "SNP.dendrogram" in incl_sections
#     ) and not "SNP.dendrogram" in excl_sections
#     if do_snp_clustering:
#         extra_sample_def = config_extract(
#             ("SNP_comparison", "extra_samples"),
#             report_settings,
#             config["reports"]["_default_"],
#         )
#         ids = collect_SNP_cluster_ids(sample_id, extra_sample_def, sample_data_full)
#         # VCF has both filtered & unfiltered SNPs now
#         # if not config_extract(('SNP_comparison', 'ignore_filter'), report_settings, config['reports']['_default_']):
#         sample_files += expand(
#             [
#                 os.path.join(
#                     DATAPATH, "{ids}", "{ids}.annotated-SNP-data.{filter}-filter.vcf.gz"
#                 )
#             ],
#             ids=ids,
#             filter=get_tool_filter_settings(
#                 f"report:{wildcards.report}:SNP_comparison"
#             ),
#         )
# 
#     if wildcards.ext == "pdf":
#         sample_files += [os.path.join(LOGPATH, "report", "_latex_installation_check")]
# 
#     return sample_files
# 

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

# rule html_report:
#     input: 
#         snp_vcf = cnv_vcf_input_function("settings:CNV_processing:call_processing"),
#         cnv_vcf = os.path.join(DATAPATH, "{sample_id}", "{sample_id}.combined-cnv-calls.vcf.gz"),
#         summary_xlsx = os.path.join(DATAPATH, "{sample_id}", "{sample_id}.summary-stats.xlsx"),
#         ref_snp_vcf = get_ref_input_function(cnv_vcf_input_function("settings:CNV_processing:call_processing")),
#         ref_cnv_vcf = get_ref_input_function('combined-cnv-calls.vcf.gz'),
#     output: 
#         os.path.join(DATAPATH, "{sample_id}", "{sample_id}.{report_name}.html"),
#         # plots=directory(os.path.join(DATAPATH,"{sample_id}","{report_name}_images")),
#     threads: get_tool_resource("knitr", "threads")
#     params:
#         # config = get_report_config
#         out_format = 'html',
#         version = __version__
#     resources:
#         runtime=get_tool_resource("knitr", "runtime"),
#         mem_mb=get_tool_resource("knitr", "memory"),
#         partition=get_tool_resource("knitr", "partition"),
#     script:
#         "../scripts/report_template.Rmd"


def get_config_delta(wildcards, compare_on=('evaluation_settings', 'settings', 'report_settings')):
    default_config = load_config(
        importlib.resources.files(STEM_CNV_CHECK).joinpath('control_files','default_config.yaml'),
        defaults=False
    )
    # Check if report name exists in default!
    default_config['report_settings'] = deep_update(
        default_config['reports'][wildcards.report] if wildcards.report in default_config['reports'] else {},
        default_config['reports']['_default_']
        
    )
    config['report_settings'] = deep_update(
        config['reports'][wildcards.report],
        config['reports']['_default_']
        
    )
    ddiff = DeepDiff(default_config, config, include_paths=compare_on, verbose_level=2)
    return ddiff.to_dict()
   

rule knit_report:
    input:
        snp_vcf = cnv_vcf_input_function("settings:CNV_processing:call_processing"),
        cnv_vcf = os.path.join(DATAPATH, "{sample_id}", "{sample_id}.combined-cnv-calls.vcf.gz"),
        summary_xlsx = os.path.join(DATAPATH, "{sample_id}", "{sample_id}.summary-stats.xlsx"),
        ref_snp_vcf = get_ref_input_function(cnv_vcf_input_function("settings:CNV_processing:call_processing")),
        ref_cnv_vcf = get_ref_input_function('combined-cnv-calls.vcf.gz'),
    output:
        report=os.path.join(DATAPATH, "{sample_id}", "{sample_id}.{report}.{ext}"),
        plots=directory(os.path.join(DATAPATH, "{sample_id}", "{report}-{ext}_images"))
    wildcard_constraints:
        ext="(pdf|html)",
    threads: get_tool_resource("knitr", "threads")
    params:
        report_config = lambda wildcards: deep_update(config['reports'][wildcards.report], config['reports']['_default_']),
        version = __version__,
        config_delta = get_config_delta,
    resources:
        runtime=get_tool_resource("knitr", "runtime"),
        mem_mb=get_tool_resource("knitr", "memory"),
        partition=get_tool_resource("knitr", "partition"),
    log:
        err=os.path.join(LOGPATH, "report", "{sample_id}", "{report}-{ext}.error.log"),
        out=os.path.join(LOGPATH, "report", "{sample_id}", "{report}-{ext}.out.log"),
    conda:
        "../envs/general-R.yaml"
    script:
        "../scripts/knit_report.R"
