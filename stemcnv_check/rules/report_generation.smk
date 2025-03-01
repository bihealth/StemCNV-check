import importlib.resources
import os
import ruamel.yaml as ruamel_yaml
from loguru import logger as logging
from deepdiff import DeepDiff
from pydantic.v1.utils import deep_update
from stemcnv_check.helpers import load_config, collect_SNP_cluster_ids, get_global_file
from stemcnv_check import __version__, STEM_CNV_CHECK

def get_penncnv_log_input(wildcards):
    sample_id = get_sample_info(wildcards)['Sample_ID']
    sex = get_sample_info(wildcards)['Sex']
    return expand(
        os.path.join(DATAPATH, sample_id, "extra_files", "PennCNV.{chrs}.error.log"),
        chrs=["auto", "chrx"] + (["chry"] if sex == "m" else []),
    )

def get_extra_snp_input_files(wildcards):
    clustering_config = config['settings']['SNV_analysis']['SNP_clustering']
    sample_id = get_sample_info(wildcards)['Sample_ID']
    # ref_id = get_sample_info(wildcards)['Reference_Sample']

    ids = set(collect_SNP_cluster_ids(sample_id, clustering_config, sample_data_df))
    # if ref_id:
    #     ids.add(ref_id)

    return expand(
        [
            os.path.join(
                DATAPATH,"{ids}","{ids}.annotated-SNP-data.{filter}-filter.vcf.gz"
            )
        ],
        ids=ids,
        filter=get_tool_filter_settings("SNV_analysis")
    )

rule make_summary_table:
    input:
        gencall_stats = os.path.join(DATAPATH, "{sample_id}", "extra_files", "{sample_id}.gencall-stats.txt"), 
        cnv_vcf = os.path.join(DATAPATH, "{sample_id}", "{sample_id}.CNV_calls.combined-annotated.vcf.gz"),
        penncnv_logs = get_penncnv_log_input,
        snv_analysis = os.path.join(DATAPATH,"{sample_id}","{sample_id}.SNV-analysis.xlsx"),
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
        err=os.path.join(LOGPATH, "summary_stats", "{sample_id}", "summary-stats.error.log"),
    conda:
        "../envs/general-R.yaml"
    script:
        "../scripts/summarise_stats.R"


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


def get_config_delta(wildcards, compare_on=('evaluation_settings', 'settings', 'report_settings')):
    default_config_path = importlib.resources.files(STEM_CNV_CHECK).joinpath('control_files','default_config.yaml')
    yaml = ruamel_yaml.YAML(typ='safe')
    with open(default_config_path) as f:
        default_config = yaml.load(f)
    # Check if report name exists in default!
    default_config['report_settings'] = deep_update(
        default_config['reports']['_default_'],
        default_config['reports'][wildcards.report] if wildcards.report in default_config['reports'] else {}
    )
    config['report_settings'] = deep_update(
        config['reports']['_default_'],
        config['reports'][wildcards.report]
    )
    ddiff = DeepDiff(default_config, config, include_paths=compare_on, verbose_level=2)
    return ddiff.to_dict()
   

rule knit_report:
    input:
        snp_vcf = snp_vcf_input_function("settings:CNV_processing:call_processing"),
        cnv_vcf = os.path.join(DATAPATH, "{sample_id}", "{sample_id}.CNV_calls.combined-annotated.vcf.gz"),
        summary_xlsx = os.path.join(DATAPATH, "{sample_id}", "{sample_id}.summary-stats.xlsx"),
        snv_analysis = os.path.join(DATAPATH,"{sample_id}","{sample_id}.SNV-analysis.xlsx"),
        ref_snp_vcf = get_ref_input_function(
            f"annotated-SNP-data.{get_tool_filter_settings('settings:CNV_processing:call_processing')}-filter.vcf.gz"
            ),
        ref_cnv_vcf = get_ref_input_function('CNV_calls.combined-annotated.vcf.gz'),
    output:
        report=os.path.join(DATAPATH, "{sample_id}", "{sample_id}.{report}.{ext}"),
        plots=directory(os.path.join(DATAPATH, "{sample_id}", "{report}-{ext}_images"))
    wildcard_constraints:
        ext="(pdf|html)",
    threads: get_tool_resource("knitr", "threads")
    params:
        report_config = lambda wildcards: deep_update(config['reports']['_default_'], config['reports'][wildcards.report]),
        version = __version__,
        config_delta = get_config_delta,
        gtf_file=lambda wildcards: get_global_file(
            'gtf', get_static_input('genome_version')(wildcards), config['global_settings'], config['cache_path']
        ),
        ginfo_file=lambda wildcards: get_global_file(
            'genome_info', get_static_input('genome_version')(wildcards), config['global_settings'], config['cache_path']
        ),
        dosage_file=lambda wildcards: get_global_file(
            'dosage_scores', get_static_input('genome_version')(wildcards), config['global_settings'], config['cache_path']
        ),
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
