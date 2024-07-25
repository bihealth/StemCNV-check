import os
import importlib.resources
from stemcnv_check import STEM_CNV_CHECK

rule filter_snp_vcf:
  input: os.path.join(DATAPATH, "{sample_id}", "{sample_id}.unprocessed.vcf")
  output: os.path.join(DATAPATH, "{sample_id}", "{sample_id}.processed-SNP-data.{filter}-filter.vcf")
  threads: get_tool_resource('filter_snp_vcf', 'threads')
  resources:
    runtime=get_tool_resource('filter_snp_vcf', 'runtime'),
    mem_mb=get_tool_resource('filter_snp_vcf', 'memory'),
    partition=get_tool_resource('filter_snp_vcf', 'partition')
  log:
    err=os.path.join(LOGPATH, "filter_snp_vcf", "{sample_id}", "{filter}.error.log"),
    #out=os.path.join(LOGPATH, "filter_snp_vcf", "{sample_id}", "{filter}.out.log")
  # conda:
  #   importlib.resources.files(STEM_CNV_CHECK).joinpath("envs","python-vcf.yaml")
  script:
    '../scripts/filter_snp_vcf.py'
    
  