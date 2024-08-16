import os

rule prep_PennCNV_sexfile:
    input: cnv_vcf_input_function('PennCNV')
    output: temp(os.path.join(DATAPATH,"{sample_id}","{sample_id}.penncnv.sexfile.txt"))
    log: os.path.join(LOGPATH,"PennCNV","{sample_id}","sexfile.log")
    params: 
        sex_info = lambda wildcards: get_ref_id(wildcards, True),
        sample_docker_path = lambda wildcards: fix_container_path(os.path.join(DATAPATH, wildcards.sample_id, f"{wildcards.sample_id}.penncnv.input.tsv"),'data'),
    shell:
        #TODO check if clause
        'echo -e "{params.sample_docker_path}\t{params.sex_info[2]}" > {output}; '
        # This can add reference sex, which is however not needed at the moment
        # 'if [[ "{sex_info[1]}" != "" ]]; then echo -e "{DATAPATH}/{sex_info[1]}/{sex_info[1]}.penncnv.input.tsv}\t{sex_info[3]}" >> {output}; fi'


# - extract tsv SNP file from vcf
rule prep_PennCNV_input:
    input:
        vcf=cnv_vcf_input_function('PennCNV')
    output:
        tsv=temp(os.path.join(DATAPATH,"{sample_id}","{sample_id}.penncnv.input.tsv"))
    log: 
        os.path.join(LOGPATH,"PennCNV","{sample_id}","input.log")
    conda:
        importlib.resources.files(STEM_CNV_CHECK).joinpath("envs","vembrane.yaml")
    params:
        filter = get_tool_filter_settings('PennCNV')
    #TODO: need to remove pseudo-autosomal regions form X & Y
    shell:
        'vembrane filter \'"PASS" in FILTER\' {input.vcf} 2> {log}|'
        'vembrane table --header \'Name, Chr, Position, B Allele Freq, Log R Ratio\''
        ' --long \'ID, CHROM, POS, FORMAT["BAF"][SAMPLE], FORMAT["LRR"][SAMPLE]\' > {output.tsv} 2>> {log}'
# # Alternative:
#     wrapper:
#         "v3.14.1/bio/vembrane/table"

# - run auto, x & Y calling
rule run_PennCNV:
    input:
        tsv=os.path.join(DATAPATH,"{sample_id}","{sample_id}.penncnv.input.tsv"),
        sexfile=os.path.join(DATAPATH,"{sample_id}","{sample_id}.penncnv.sexfile.txt"),
        pfb=config['static_data']['penncnv_pfb_file'],
        gcmodel=config['static_data']['penncnv_GCmodel_file']
    output:
        tsv=temp(os.path.join(DATAPATH,"{sample_id}","{sample_id}.penncnv-{chr}.tsv")),
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
        #LOH calling is optional, also disabled for male sex chromosomes
        do_loh=lambda wildcards: '' if (wildcards.chr != 'auto' and get_ref_id(wildcards,True)[2] == 'm' and 
                                        config['settings']['PennCNV']['enable_LOH_calls']) else '-loh',
        snakedir = fix_container_path(SNAKEDIR, 'snakedir'),
        tsvout=lambda wildcards: fix_container_path(os.path.join(DATAPATH, wildcards.sample_id, wildcards.sample_id+".penncnv-"+wildcards.chr+".tsv"),'data'),
        tsvin=lambda wildcards: fix_container_path(os.path.join(DATAPATH, wildcards.sample_id, wildcards.sample_id+".penncnv.input.tsv"),'data'),
        pfb=fix_container_path(config['static_data']['penncnv_pfb_file'],'static'),
        sexfile=lambda wildcards: '' if wildcards.chr == 'auto' else '--sexfile ' + str(fix_container_path(os.path.join(DATAPATH, wildcards.sample_id, wildcards.sample_id+".penncnv.sexfile.txt"),'data')),
        gcmodel=fix_container_path(config['static_data']['penncnv_GCmodel_file'],'static'),
        logerr=lambda wildcards:fix_container_path(os.path.join(LOGPATH,"PennCNV", wildcards.sample_id, wildcards.chr+".error.log"),'logs'),
        logout=lambda wildcards:fix_container_path(os.path.join(LOGPATH,"PennCNV", wildcards.sample_id, wildcards.chr+".out.log"),'logs')
    log:
        err=os.path.join(LOGPATH,"PennCNV","{sample_id}","{chr}.error.log"),
        out=os.path.join(LOGPATH,"PennCNV","{sample_id}","{chr}.out.log")
    container:
        "docker://genomicslab/penncnv"
    shell:
        '/home/user/PennCNV/detect_cnv.pl -test {params.do_loh} {params.chrom} -confidence '
        '-hmm {params.snakedir}/supplemental-files/hhall_loh.hmm -pfb {params.pfb} -gcmodel {params.gcmodel} '
        '{params.sexfile} {params.tsvin} -out {params.tsvout} > {params.logout} 2> {params.logerr}'


def get_penncnv_output(wildcards, files = 'tsv'):
    _, _, sex, _ = get_ref_id(wildcards, True) 
    chrs = ['auto','chrx']
    if sex == 'm': chrs.append('chry')
    
    if files == 'tsv':
        return expand(
            os.path.join(DATAPATH,"{sample_id}","{sample_id}.penncnv-{chr}.tsv"),
            chr=chrs, sample_id=wildcards.sample_id
        )
    elif files == 'log':
        return expand(
            os.path.join(LOGPATH,"PennCNV","{sample_id}","{chr}.error.log"),
            chr=chrs, sample_id=wildcards.sample_id
        )
    else:
        raise ValueError("Invalid file type requested: {}".format(files))


# - combine calls and write out as vcf
rule combined_PennCNV_output:
    input:
        vcf=cnv_vcf_input_function('PennCNV'),
        tsvs=get_penncnv_output,
        #logs=get_penncnv_output(wildcards, 'log')
    output:
        vcf=os.path.join(DATAPATH,"{sample_id}","{sample_id}.CNV_calls.PennCNV.vcf.gz"),
        # TODO: how to best handle this? > separate rule ?!
        # stats=os.path.join(DATAPATH,"{sample_id}","{sample_id}.CNV_calls.penncnv.stats.tsv")
    log:
        err=os.path.join(LOGPATH,"PennCNV","{sample_id}","combine.error.log"),
    conda:
        importlib.resources.files(STEM_CNV_CHECK).joinpath("envs","general-R.yaml")
    script: 
        "../scripts/combine_penncnv.R"