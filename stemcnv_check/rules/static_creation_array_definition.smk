# -*- coding: utf-8 -*-
import importlib.resources
import os
from pathlib import Path
import tempfile
from stemcnv_check import STEM_CNV_CHECK, ENSEMBL_RELEASE, mehari_db_version

DOWNLOAD_DIR = config["TMPDIR"] if "TMPDIR" in config else tempfile.mkdtemp()
GENOME = config["genome"]
ARRAY = config["array_name"]
# ================================================================


def fix_container_path(path_in, bound_to):
    path_in = Path(path_in)

    if bound_to in ("static", ARRAY):
        rel_path = path_in.name
    else:
        local_target = {
            "snakedir": Path(importlib.resources.files(STEM_CNV_CHECK)),
            "tmp": Path(DOWNLOAD_DIR),
        }[bound_to].absolute()
        rel_path = path_in.absolute().relative_to(local_target)

    return Path("/outside/") / bound_to / rel_path


# Note: PennCNV does not seem to work with UCSC chromosome style in PFB file
rule create_pfb_from_vcf:
    input:
        config["vcf_input_file"],
    output:
        config["penncnv_pfb_file"],
    conda:
        "../envs/general-R.yaml"
    shell:
        """
Rscript - << 'EOF'
suppressMessages(library(tidyverse))
suppressMessages(library(vcfR))
suppressMessages(library(GenomeInfoDb))

snp.vcf <- read.vcfR('{input}', verbose = F)

vcf.info <- vcfR2tidy(
        snp.vcf, 
        info_only=T,
        info_fields=c('ALLELE_A', 'ALLELE_B', 'N_AA', 'N_AB', 'N_BB') # 'GC',
    ) %>%
    .$fix %>%
    mutate(
        A_freq = (2*N_AA + N_AB) / (2*(N_AA + N_AB + N_BB)),
        #This is B_allele frequency, or 'population frequency B-allele'. PennCNV needs this
        PFB = (2*N_BB + N_AB) / (2*(N_AA + N_AB + N_BB)),
        #This is AF/Alternate allele frequency. Some bcftool modules need this
        Alt_freq = ifelse(ALLELE_B == 1, PFB, A_freq)
    )

# ensure that CHROM uses NCBI/ENSEMBL styles
seqlevelsStyle(vcf.info$CHROM) <- 'NCBI'

vcf.info %>%
    mutate(Chr = factor(CHROM, levels=genomeStyles()$Homo_sapiens$NCBI)) %>%
    select(ID, Chr, POS, PFB) %>%
    dplyr::rename(Name = ID, Position = POS) %>%
    filter(if_all(everything(), ~!is.na(.))) %>%
    arrange(Chr, Position) %>%
    write_tsv('{output}')
EOF
"""


# FUTURE:
# also get the gc5base.bw file from UCSCand make the GC model from that?
# -> PennCNV comes with a wig2gc5base python script (though that has a hard coded 'source' file in it?
rule create_gcmodel_file:
    input:
        config["penncnv_pfb_file"],
    output:
        config["penncnv_GCmodel_file"],
    params:
        download_path=fix_container_path(DOWNLOAD_DIR, "tmp"),
        input_path=fix_container_path(config["penncnv_pfb_file"], ARRAY),
        output_path=fix_container_path(config["penncnv_GCmodel_file"], ARRAY)
    container:
        "docker://genomicslab/penncnv"
    shell:
        """
        if [[ -f /home/user/PennCNV/gc_file/{GENOME}.gc5Base.txt ]]; then
            ln -s /home/user/PennCNV/gc_file/{GENOME}.gc5Base.txt {params.download_path}/{GENOME}.gc5Base.txt
        else
            gunzip -c /home/user/PennCNV/gc_file/{GENOME}.gc5Base.txt.gz > {params.download_path}/{GENOME}.gc5Base.txt 
        fi
        /home/user/PennCNV/cal_gc_snp.pl {params.download_path}/{GENOME}.gc5Base.txt {params.input_path} -out {params.output_path}
        """


rule create_array_info_files:
    input:
        pfb=config["penncnv_pfb_file"],
        chromInfo=config["genomeInfo"],
    output:
        density=config["array_density_file"],
        gaps=config["array_gaps_file"],
    params:
        min_gap_size=config["min_gap_size"],
        density_windows=config["density_windows"],
        genome=config["genome"],
    conda:
        "../envs/general-R.yaml"
    shell:
        """
Rscript - << 'EOF'
suppressMessages(library(tidyverse))
suppressMessages(library(plyranges))
suppressMessages(library(GenomeInfoDb))

# The 'getChromInfoFromUCSC' function can break if UCSC does updates faster than Bioconductor
# therefore using the chromInfo.txt file from UCSC instead
chromInfo <- read_tsv('{input.chromInfo}', show_col_types = FALSE) %>%
    mutate(chr = factor(chr, levels=unique(chr))) %>%
    group_by(chr) %>% summarise(seqlength = unique(size))
ginfo <- Seqinfo(
    seqnames=chromInfo$chr %>% as.character(),
    seqlengths=chromInfo$seqlength,
    isCircular= chromInfo$chr %in% c('chrMT', 'chrM'),
    genome='{params.genome}'
)

# Read pfb to get probe positions, pfb uses NCBI chromosome style
array <- read_tsv('{input.pfb}', show_col_types = FALSE) %>%
  as_granges(seqnames = Chr, start = Position, width = 1) %>%
  #Need to filter duplicate positions, as they screw up the density
  reduce_ranges(n_probes = n())
  
seqlevelsStyle(array) <- 'UCSC'
array <- filter(array, seqnames %in% seqnames(ginfo))
seqlevels(array) <- seqlevels(ginfo)
seqinfo(array) <- ginfo

# Get gaps
gaps <- complement_ranges(array) %>%
    mutate(
        chrom_start = start == 1,
        chrom_end = end == seqlengths(.)[as.character(seqnames)],
    ) %>%
    filter(!chrom_end & !chrom_start) %>%
    mutate(
        width_mean_sd = mean(width(.))+sd(width(.)),
        gap_min_size = ifelse(
            '{params.min_gap_size}' == 'auto-array',
            unique(width_mean_sd),
            as.numeric('{params.min_gap_size}')
        )
    ) %>%
    filter(width(.) >= gap_min_size)
write_bed(gaps, '{output.gaps}')

# calculate probe density 
# by tiling the whole genome we have even windows, and mean & median are similar
density <- chromInfo %>% 
    as_granges(seqnames = chr, width = seqlength, start = 1) %>%
    tile_ranges({params.density_windows}) %>%
    mutate(score = count_overlaps(., array) / width(.))
write_bed(density, '{output.density}')

EOF
"""


