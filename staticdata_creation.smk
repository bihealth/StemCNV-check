# -*- coding: utf-8 -*-

import os
import tempfile
import yaml
from scripts.py_helpers import *
from scripts.py_exceptions import *

SNAKEDIR = config['snakedir'] if 'snakedir' in config else os.path.dirname(os.path.realpath(__file__)) #Defined by wrapper
DOWNLOAD_DIR = config['TMPDIR'] if 'TMPDIR' in config else tempfile.mkdtemp()
GENOME = config['genome']

# ================================================================

rule all:
    input: [config['chrominfo_outname'], config['pfb_output_name'], config['gcmodel_output_name']]
#FUTURE:
# a nicer naming scheme for pfb & gcmodel files would require somewhat standardised naming scheme for different arrays
# without that it's easier to have the naming be user/wrapper defined

rule create_pfb_from_vcf:
    input: config['vcf_input_file']
    output: config['pfb_outname']
    shell: "Rscript {SNAKEDIR}/scripts/make_PFB_from_vcf.R {input} {output}"

#FUTURE:
#also get the gc5base.bw file from UCSCand make the GC model from that?
# -> PennCNV comes with a wig2gc5base python script (though that has a hard coded 'source' file in it?
rule create_gcmodel_file:
    input: config['pfb_outname']
    output: config['gcmodel_outname']
    params:
        penncnv_path = '/home/user/PennCNV' if config['use_singularity'] else '$CONDA_PREFIX/pipeline/PennCNV-1.0.5'
    container:
        "docker://genomicslab/penncnv"
    shell:
        "gunzip {params.penncnv_path}/gc_file/{GENOME}.gc5Base.txt.gz; {params.penncnv_path}/cal_gc_snp.pl {params.penncnv_path}/gc_file/{GENOME}.gc5Base.txt {input} -out {output}"

#TODO not all files exist for all Genome versions!
# -> might need a check that download worked?
rule ucsc_goldenpath_download:
    output: temp("{DOWNLOAD_DIR}/{genome}.{filename}.txt")
    wildcard_constraints:
        genome = "hg19|hg38"
    shell:
        """
        wget 'ftp://hgdownload.cse.ucsc.edu/goldenPath/{wildcards.genome}/database/{wildcards.filename}.txt.gz' -O {output}.gz 2> /dev/null
        gunzip {output}.gz
        """


rule ucsc_genome_positions:
    input:
        cytobands = os.path.join(DOWNLOAD_DIR, f"{GENOME}.cytoBand.txt"),
        centromer = os.path.join(DOWNLOAD_DIR, f"{GENOME}.centromeres.txt"),
        chrominfo = os.path.join(DOWNLOAD_DIR, f"{GENOME}.chromInfo.txt")
        #Note: the chromAlias.txt file might be useful if people use strange chr-/seqnames
    output: config['chrominfo_outname']
    shell:
        """
Rscript - << 'EOF'
library(tidyverse)
message('startup')
chrominfo <- read_tsv("{input.chrominfo}", col_names = c("chr", "size", "url")) %>% select(-url) %>% filter(str_detect(chr, "_(alt|fix)"))
chrominfo <- merge(chrominfo, by = "chr", all=T,
                read_tsv("{input.centromer}", col_names = c("bin", "chr", "start", "end", "name")) %>%
                    group_by(chr) %>% summarise(centromer_start = min(start), centromer_end = max(end)))
chrominfo <- merge(chrominfo, by = "chr", all=T,
                read_tsv("{input.cytobands}", col_names = c("chr", "band_start", "band_end", "band_name", "band_staining"))) 
write_tsv(chrominfo, "{output}")
EOF
        """


