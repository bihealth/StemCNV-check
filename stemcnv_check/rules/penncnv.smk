import os


# New rules needed:
# ?? write penncnv sex file (only for sample)
# - extract tsv SNP file from vcf
# - run auto, x & Y calling
# - combine calls and write out as vcf

if config['use_singularity']:
  rule run_PennCNV:
    input:
      tsv=os.path.join(DATAPATH,"{sample_id}","{{sample_id}}.filtered-data-{0}.tsv".format(get_tool_filter_settings('PennCNV'))),
      pfb=config['static_data']['penncnv_pfb_file'],
      gcmodel=config['static_data']['penncnv_GCmodel_file']
    output:
      tsv=os.path.join(DATAPATH,"{sample_id}","{sample_id}.penncnv-{chr}.tsv"),
      err=os.path.join(LOGPATH,"PennCNV","{sample_id}","{chr}.error.log")
    threads: get_tool_resource('PennCNV','threads')
    resources:
      runtime=get_tool_resource('PennCNV','runtime'),
      mem_mb=get_tool_resource('PennCNV','memory'),
      partition=get_tool_resource('PennCNV','partition')
    wildcard_constraints:
      chr='chrx|chry|auto'
    params:
      filter=get_tool_filter_settings('PennCNV'),
      chrom=lambda wildcards: '' if wildcards.chr == 'auto' else '-' + wildcards.chr,
      #Male sex chromosomes can't have LOH, but PennCNV does not exclude it if run with -loh
      do_loh=lambda wildcards: '' if (wildcards.chr != 'auto' and get_ref_id(wildcards,True)[2] == 'm') else '-loh',
      snakedir = fix_container_path(SNAKEDIR, 'snakedir'),
      tsvout=lambda wildcards: fix_container_path(os.path.join(DATAPATH, wildcards.sample_id, wildcards.sample_id+".penncnv-"+wildcards.chr+".tsv"),'data'),
      tsvin=lambda wildcards: fix_container_path(os.path.join(DATAPATH, wildcards.sample_id, wildcards.sample_id+".filtered-data-{}.tsv".format(get_tool_filter_settings('PennCNV'))),'data'),
      pfb=fix_container_path(config['static_data']['penncnv_pfb_file'],'static'),
      gcmodel=fix_container_path(config['static_data']['penncnv_GCmodel_file'],'static'),
      logerr=lambda wildcards:fix_container_path(os.path.join(LOGPATH,"PennCNV", wildcards.sample_id, wildcards.chr+".error.log"),'logs'),
      logout=lambda wildcards:fix_container_path(os.path.join(LOGPATH,"PennCNV", wildcards.sample_id, wildcards.chr+".out.log"),'logs')
    log:
      err=os.path.join(LOGPATH,"PennCNV","{sample_id}","{chr}.error.log"),
      out=os.path.join(LOGPATH,"PennCNV","{sample_id}","{chr}.out.log")
    container:
      "docker://genomicslab/penncnv"
    shell:
      '/home/user/PennCNV/detect_cnv.pl -test {params.do_loh} {params.chrom} -confidence -hmm {params.snakedir}/supplemental-files/hhall_loh.hmm -pfb {params.pfb} -gcmodel {params.gcmodel} {params.tsvin} -out {params.tsvout} > {params.logout} 2> {params.logerr}'
else:
  rule run_PennCNV:
    input:
      tsv = os.path.join(DATAPATH, "{sample_id}", "{{sample_id}}.filtered-data-{0}.tsv".format(get_tool_filter_settings('PennCNV'))),
      pfb = config['static_data']['penncnv_pfb_file'],
      gcmodel = config['static_data']['penncnv_GCmodel_file']
    output:
      tsv = os.path.join(DATAPATH, "{sample_id}", "{sample_id}.penncnv-{chr}.tsv"),
      err=os.path.join(LOGPATH,"PennCNV","{sample_id}","{chr}.error.log"),
    threads: get_tool_resource('PennCNV', 'threads')
    resources:
      runtime=get_tool_resource('PennCNV', 'runtime'),
      mem_mb=get_tool_resource('PennCNV', 'memory'),
      partition=get_tool_resource('PennCNV', 'partition')
    wildcard_constraints:
      chr='chrx|chry|auto'
    params:
      filter = get_tool_filter_settings('PennCNV'),
      chrom = lambda wildcards: '' if wildcards.chr == 'auto' else '-'+wildcards.chr,
      #Male sex chromosomes can't have LOH, but PennCNV does not exclude it if run with -loh
      do_loh = lambda wildcards: '' if (wildcards.chr != 'auto' and get_ref_id(wildcards, True)[2] == 'm') else '-loh'
    log:
      err=os.path.join(LOGPATH, "PennCNV", "{sample_id}", "{chr}.error.log"),
      out=os.path.join(LOGPATH, "PennCNV", "{sample_id}", "{chr}.out.log")
    shell:
      'PennCNV_detect -test {params.do_loh} {params.chrom} -confidence -hmm {SNAKEDIR}/supplemental-files/hhall_loh.hmm -pfb {input.pfb} -gcmodel {input.gcmodel} {input.tsv} -out {output.tsv} > {log.out} 2> {log.err}'
