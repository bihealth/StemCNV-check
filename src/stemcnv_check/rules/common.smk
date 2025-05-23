import importlib.resources
import os
from pathlib import Path
import ruamel.yaml as ruamel_yaml
from stemcnv_check import STEM_CNV_CHECK
from stemcnv_check.helpers import config_extract, get_global_file, get_array_file
from stemcnv_check.exceptions import SampleConstraintError
from loguru import logger as logging

# Determine if we can pipe initial vcf files
# Do not pipe if:
# 1) running on WSL and writing to /mnt/[any windows drive]
#FIXME:
# 2) <3 cores or not enough memory for all 3 rules 
# > requires to determine available cores and set memory
if (
    config['is_wsl'] and Path(config['data_path']).absolute().match("/mnt/*") 
):
    pipe_or_temp_function = temp
    logging.warning(
        "Piping is disabled because the pipeline is running on WSL and writing to a Windows drive, "
        "this may result in slower performance. "
    )
else:
    pipe_or_temp_function = pipe

def get_sample_info(wildcards):
    info = sample_data_df.loc[wildcards.sample_id].to_dict()
    # ensure Sex is encoded as 'm' or 'f'
    info['Sex'] = info['Sex'][0].lower()
    return info

def fix_container_path(path_in, bound_to):
    path_in = Path(path_in)

    local_targets = {
        "data": Path(DATAPATH),
        "rawdata": Path(IDAT_INPUT),
        "logs": Path(LOGPATH),
        "snakedir": Path(importlib.resources.files(STEM_CNV_CHECK)),
    }
    # Get relative paths for local folders of /created by the pipeline 
    if bound_to in local_targets.keys():
        local_target = local_targets[bound_to].absolute()
        rel_path = path_in.absolute().relative_to(local_target)
    # static data files could be distributed over different systems, so are bound 1-by-1
    # Need to allow different 'bound_to' values for different array static files
    else:
        rel_path = path_in.name

    return Path("/outside/") / bound_to / rel_path


def get_tool_filter_settings(tool):
    if tool.split(":")[0] == "report":
        report_settings = config["reports"][tool.split(":")[1]]
        out = config_extract(
            (tool.split(":")[2], "probe_filter_settings"),
            report_settings,
            config["reports"]["_default_"],
        )
    elif tool.count(":") == 2 and tool.split(":")[1] == "CNV_processing":
        out = config["settings"]["CNV_processing"]["call_processing"]["probe_filter_settings"]
    else:
        out = config["settings"][tool]["probe_filter_settings"]
        if tool == 'SNV_analysis':
            out = out if out != "none" else '_default_'
    if out == "_default_":
        out = config["settings"]["default_probe_filter_set"]
    return out


def get_tool_resource(tool, resource):
    if not resource in ("threads", "memory", "runtime", "partition", "cmd-line-params"):
        raise KeyError(f"This resource can not be defined: {resource}")
    if not tool in config["tools"] or not resource in config["tools"][tool]:
        return config["tools"]["_default_"][resource]
    else:
        return config["tools"][tool][resource]


def get_genome_fasta(wildcards):
    """Get the correct genome fasta file"""

    chip = get_sample_info(wildcards.sample_id)['Array_Name']
    genome_version = config['array_definition'][chip]['genome_version']
    
    return get_global_file('fasta', genome_version, config['global_settings'], config['cache_path'])

def get_static_input(type):
    """Get the correct static files matching array and genome version"""
    
    if type not in (
        'fasta', 'gtf', 'genome_info', 'mehari_txdb',
        'genome_version', 'bpm_manifest_file', 'egt_cluster_file', 'csv_manifest_file',
        'penncnv_pfb_file', 'penncnv_GCmodel_file', 'array_density_file', 'array_gaps_file'
    ):
        raise ValueError(f"Invalid filetype: {type}")
    
    def input_function(wildcards):

        if hasattr(wildcards, 'sample_id'):
            array = get_sample_info(wildcards)['Array_Name']
        elif hasattr(wildcards, 'sentrix_name'):
            array = set(sample_data_df.loc[sample_data_df['Chip_Name'] == wildcards.sentrix_name]['Array_Name'].values)
            if len(array) != 1:
                raise ValueError(f"Multiple array types found for sentrix_name {wildcards.sentrix_name}: {array}")
            array = array.pop()
        else:
            raise ValueError("No sample_id or other matching attribute found in wildcards")

        genome = config['array_definition'][array]['genome_version']
    
        res_dict = {
            file_type: get_global_file(file_type, genome, config['global_settings'], config['cache_path'])
            for file_type in ('fasta', 'gtf', 'genome_info', 'mehari_txdb')
        }
        res_dict.update({
            k: v for k, v in config['array_definition'][array].items()
            if not k.startswith('array') and not k.startswith('penncnv')
        })
        res_dict.update({
            key: get_array_file(key, array, config, config['cache_path'])
            for key in config['array_definition'][array].keys()
            if key.startswith('array') or key.startswith('penncnv')
        })
        
        return res_dict[type]
    
    return input_function


def snp_vcf_input_function(tool):
    return lambda wildcards: os.path.join(
        DATAPATH,
        wildcards.sample_id,
        # f"{wildcards.sample_id}.processed-SNP-data.{get_tool_filter_settings(tool)}-filter.vcf",
        f"{wildcards.sample_id}.annotated-SNP-data.{get_tool_filter_settings(tool)}-filter.vcf.gz"
    )

def get_ref_input_function(input_file_type):
    # Check if file_pattern is already a (input) function
    if callable(input_file_type):
        get_file_pattern = input_file_type
    else:
        get_file_pattern = lambda wildcards: os.path.join(DATAPATH, "{sample_id}", f"{{sample_id}}.{input_file_type}")
    def input_function(wildcards):
        sample_id = get_sample_info(wildcards)['Sample_ID']
        ref_id = get_sample_info(wildcards)['Reference_Sample']
        if ref_id:
            return expand(get_file_pattern(wildcards), sample_id = ref_id)
        else:
            return []
    return input_function

def fix_call_label_order(call_labels):
    """Fix the order of the call labels, so that defaults from label_name_definitions.yaml file come last"""
    yaml = ruamel_yaml.YAML(typ='safe')
    with open(importlib.resources.files(STEM_CNV_CHECK).joinpath('control_files','label_name_definitions.yaml')) as f:
        defined_call_labels = yaml.load(f)['CNV_labels']

    labels_ordered = {k: v for k, v in call_labels.items() if k not in defined_call_labels}
    labels_ordered.update(call_labels)

    return labels_ordered
