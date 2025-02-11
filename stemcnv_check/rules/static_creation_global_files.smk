# -*- coding: utf-8 -*-
import os
from loguru import logger as logging
from stemcnv_check import ENSEMBL_RELEASE
from stemcnv_check.helpers import get_global_file

GENOME = config["genome"]
DOWNLOAD_DIR = config["TMPDIR"] if "TMPDIR" in config else tempfile.mkdtemp()

rule download_gencode_gtf_v45:
    output:
        get_global_file('gtf', GENOME, config['global_settings'], config['cache_path'], False)
    params:
        ftp_base="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/",
        release_path=(
            "gencode.v45.basic.annotation.gtf.gz"
            if GENOME == "hg38"
            else "GRCh37_mapping/gencode.v45lift37.basic.annotation.gtf.gz"
        ),
    shell:
        "wget {params.ftp_base}/{params.release_path} -O {output} 2> /dev/null"

rule download_ensmble_fasta:
    output:
        get_global_file('fasta', GENOME, config['global_settings'], config['cache_path'], False).removesuffix(".gz")
    params:
        species="homo_sapiens",
        datatype="dna",
        release=ENSEMBL_RELEASE,
        build='GRCh37' if GENOME in ('hg19', 'GRCh37') else 'GRCh38',
    cache: "omit-software"  # save space and time with between workflow caching (see docs)
    wrapper:
        "v5.5.0/bio/reference/ensembl-sequence"

rule bgzip_fasta:
    input:
        get_global_file('fasta',GENOME,config['global_settings'],config['cache_path'],False).removesuffix(".gz")
    output:
        get_global_file('fasta', GENOME, config['global_settings'], config['cache_path'], False)
    wrapper:
        "v5.5.0/bio/bgzip"


rule download_mehari_ensembl_db:
    output:
        get_global_file('mehari_txdb', GENOME, config['global_settings'],config['cache_path'], False)
    wildcard_constraints: 
        genome = 'GRCh37|GRCh38',
        MEHARI_DB_VERSION = '[0-9]\\.[0-9]\\.[0-9]'
    shell:
        "wget https://github.com/varfish-org/mehari-data-tx/releases/download/v{wildcards.MEHARI_DB_VERSION}/$(basename {output}) -O {output}"


rule ucsc_goldenpath_download:
    output:
        temp("{DOWNLOAD_DIR}/{genome}.{filename}.txt"),
    wildcard_constraints:
        genome="hg19|hg38",
    shell:
        """
        wget 'https://hgdownload.cse.ucsc.edu/goldenPath/{wildcards.genome}/database/{wildcards.filename}.txt.gz' -O {output}.gz 2> /dev/null
        gunzip {output}.gz
        """


rule create_genome_info_file:
    input:
        cytobands=ancient(os.path.join(DOWNLOAD_DIR, f"{GENOME}.cytoBand.txt")),
        #cytoBandIdeo.txt has info for non-assmebled chromosomes
        #centromer = os.path.join(DOWNLOAD_DIR, f"{GENOME}.centromeres.txt"), #Only exists for hg38
        chrominfo=ancient(os.path.join(DOWNLOAD_DIR, f"{GENOME}.chromInfo.txt")),
        #Note: the chromAlias.txt file might be useful if people use strange chr-/seqnames
    output:
        get_global_file('genome_info', GENOME, config['global_settings'], config['cache_path'], False)
        # config["genomeInfo"],
    conda:
        "../envs/general-R.yaml"
    shell:
        """
Rscript - << 'EOF'
suppressMessages(library(tidyverse))
suppressMessages(library(GenomeInfoDb))

style <- genomeStyles()$Homo_sapiens$UCSC %>% str_replace('chrM', 'chrMT')

chrominfo <- read_tsv("{input.chrominfo}", col_names = c("chr", "size", "url"), show_col_types = FALSE) %>% 
    select(-url) %>% 
    # filter(chr != 'chrMT') %>%
    # mutate(chr = str_replace(chr, 'chrM$', 'chrMT')) %>%
    filter(chr %in% style)
# The 'acen' labelled cytobands are the centromere areas (but much larger than what the centromere file gives)
chrominfo <- left_join(
    chrominfo,
    read_tsv(
        "{input.cytobands}",
        show_col_types = FALSE, 
        col_names = c("chr", "band_start", "band_end", "band_name", "band_staining")
    ),
    by = "chr"
) %>%
    mutate(
        centromer = band_staining == "acen",
        chr = factor(chr, style)
    ) %>%
    arrange(chr, band_start) %>%
    # The hg19 cytoband file does not include chrMT, so no band_start there
    filter(!is.na(band_start) | chr == 'chrM')
write_tsv(chrominfo, "{output}")
EOF
        """

rule download_dosage_sensivity_data:
    output: 
        get_global_file('dosage_scores', GENOME, config['global_settings'], config['cache_path'], False)
    shell:
        """
        wget 'https://zenodo.org/records/6347673/files/Collins_rCNV_2022.dosage_sensitivity_scores.tsv.gz' -O {output} 2> /dev/null
        """