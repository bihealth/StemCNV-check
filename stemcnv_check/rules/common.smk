# Makestatic Data :
# def fix_container_path(path_in, bound_to):
# 
#     path_in = Path(path_in)
#     if bound_to == 'static':
#         rel_path = path_in.name
#     else:
#         local_target = {
#             'snakedir': Path(importlib.resources.files(STEM_CNV_CHECK)),
#             'tmp': Path(DOWNLOAD_DIR)
#         }[bound_to].absolute()
#         rel_path = path_in.absolute().relative_to(local_target)
# 
#     return Path('/outside/') / bound_to / rel_path

def fix_container_path(path_in, bound_to):
  path_in = Path(path_in)

  if bound_to == 'static':
    rel_path = path_in.name
  else:
    local_target = {
      'data': Path(DATAPATH),
      'rawdata': Path(IDAT_INPUT),
      'logs': Path(LOGPATH),
      'snakedir': Path(importlib.resources.files(STEM_CNV_CHECK)),
    }[bound_to].absolute()
    rel_path = path_in.absolute().relative_to(local_target)

  return Path('/outside/') / bound_to / rel_path

def get_tool_filter_settings(tool):
  if tool.split(':')[0] == 'report':
    report_settings = config['reports'][tool.split(':')[1]]
    out = config_extract((tool.split(':')[2], 'filter-settings'), report_settings, config['reports']['__default__'])
  elif tool.count(':') == 2 and tool.split(':')[1] == 'CNV_processing':
    out = config['settings']['CNV_processing']['call_processing']['filter-settings']
  else:
    out = config['settings'][tool]['filter-settings']
  if out == '__default__':
    out = config['settings']['default-filter-set']
  return out


def get_tool_resource(tool ,resource):
  if not resource in ('threads', 'memory', 'runtime', 'partition', 'cmd-line-params'):
    raise KeyError(f"This resource can not be defined: {resource}")
  if not tool in config['tools']:
    return config['tools']['__default__'][resource]
  else:
    if resource in config['tools'][tool]:
      return config['tools'][tool][resource]
    else:
      return config['tools']['__default__'][resource]

def get_genome_fasta(wildcards):
  """Get the correct genome fasta file"""
  # #FIXME: future    
  # chip = get_sample_info(wildcards.sample_id)['array_name']
  # genome = config['array_definitions'][chip]['genome_version']
  if config['genome_version'] in ('hg38', 'GRCh38'):
    return config['global_settings']['hg38_genome_fasta']
  else:
    return config['global_settings']['hg19_genome_fasta']