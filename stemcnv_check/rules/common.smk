import importlib.resources
import os
from pathlib import Path
from stemcnv_check import STEM_CNV_CHECK
from stemcnv_check.helpers import config_extract
from stemcnv_check.exceptions import SampleConstraintError

# Makestatic Data :
# def fix_container_path(path_in, bound_to):
#
#         path_in = Path(path_in)
#         if bound_to == 'static':
#                 rel_path = path_in.name
#         else:
#                 local_target = {
#                         'snakedir': Path(importlib.resources.files(STEM_CNV_CHECK)),
#                         'tmp': Path(DOWNLOAD_DIR)
#                 }[bound_to].absolute()
#                 rel_path = path_in.absolute().relative_to(local_target)
#
#         return Path('/outside/') / bound_to / rel_path


def get_ref_id(wildcards, get_sex=False):
    sample_id = wildcards.sample_id
    sex, ref_id = [(s, r) for sid, _, _, s, r in sample_data if sid == sample_id][0]
    sex = sex[0].lower()

    if ref_id:
        try:
            # Assume existing match -> wrapper should have done a check
            ref_sex = [s for sid, _, _, s, _ in sample_data if sid == ref_id][0]
            ref_sex = ref_sex[0].lower()
        except IndexError:
            # Somehow no match
            raise SampleConstraintError(
                f"Listed reference sample can not be found in sample-table: '{ref_id}'"
            )
    else:
        ref_id = False
        ref_sex = False

    if get_sex:
        return sample_id, ref_id, sex, ref_sex
    else:
        return sample_id, ref_id


def fix_container_path(path_in, bound_to):
    path_in = Path(path_in)

    if bound_to == "static":
        rel_path = path_in.name
    else:
        local_target = {
            "data": Path(DATAPATH),
            "rawdata": Path(IDAT_INPUT),
            "logs": Path(LOGPATH),
            "snakedir": Path(importlib.resources.files(STEM_CNV_CHECK)),
        }[bound_to].absolute()
        rel_path = path_in.absolute().relative_to(local_target)

    return Path("/outside/") / bound_to / rel_path


def get_tool_filter_settings(tool):
    if tool.split(":")[0] == "report":
        report_settings = config["reports"][tool.split(":")[1]]
        out = config_extract(
            (tool.split(":")[2], "filter-settings"),
            report_settings,
            config["reports"]["__default__"],
        )
    elif tool.count(":") == 2 and tool.split(":")[1] == "CNV_processing":
        out = config["settings"]["CNV_processing"]["call_processing"]["filter-settings"]
    else:
        out = config["settings"][tool]["filter-settings"]
    if out == "__default__":
        out = config["settings"]["default-filter-set"]
    return out


def get_tool_resource(tool, resource):
    if not resource in ("threads", "memory", "runtime", "partition", "cmd-line-params"):
        raise KeyError(f"This resource can not be defined: {resource}")
    if not tool in config["tools"] or not resource in config["tools"][tool]:
        return config["tools"]["__default__"][resource]
    else:
        return config["tools"][tool][resource]


def get_genome_fasta(wildcards):
    """Get the correct genome fasta file"""
    # #FIXME: future
    # chip = get_sample_info(wildcards.sample_id)['array_name']
    # genome = config['array_definitions'][chip]['genome_version']
    if config["genome_version"] in ("hg38", "GRCh38"):
        return config["global_settings"]["hg38_genome_fasta"]
    else:
        return config["global_settings"]["hg19_genome_fasta"]


def cnv_vcf_input_function(tool):
    return lambda wildcards: os.path.join(
        DATAPATH,
        wildcards.sample_id,
        f"{wildcards.sample_id}.processed-SNP-data.{get_tool_filter_settings(tool)}-filter.vcf",
        # VEP still has issues; skip for now
        # f"{wildcards.sample_id}.annotated-SNP-data.{get_tool_filter_settings(tool)}-filter.vcf.gz"
    )
