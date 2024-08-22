# -*- coding: utf-8 -*-
import importlib.resources
import os
from pathlib import Path
import tempfile
from stemcnv_check import STEM_CNV_CHECK

DOWNLOAD_DIR = config["TMPDIR"] if "TMPDIR" in config else tempfile.mkdtemp()
GENOME = config["genome"]

# ================================================================


def fix_container_path(path_in, bound_to):
    path_in = Path(path_in)

    if bound_to == "static":
        rel_path = path_in.name
    else:
        local_target = {
            "snakedir": Path(importlib.resources.files(STEM_CNV_CHECK)),
            "tmp": Path(DOWNLOAD_DIR),
        }[bound_to].absolute()
        rel_path = path_in.absolute().relative_to(local_target)

    return Path("/outside/") / bound_to / rel_path


#
# rule all:
#     input:
#         config['genomeInfo_file'],
#         config['array_gaps_file'],
#         config['array_density_file'],
#         config['penncnv_pfb_file'],
#         config['penncnv_GCmodel_file'],
#         config['genome_gtf_file']


# Note: PennCNV does not seem to work with UCSC chromosome style in PFB file
rule create_pfb_from_vcf:
    input:
        config["vcf_input_file"],
    output:
        config["penncnv_pfb_file"],
    conda:
        "../envs/general-R.yaml"
    # shell: "Rscript {SNAKEDIR}/scripts/make_PFB_from_vcf.R {input} {output}"
    shell:
        """
Rscript - << 'EOF'
suppressMessages(library(tidyverse))
suppressMessages(library(vcfR))
suppressMessages(library(GenomeInfoDb))

snp.vcf <- read.vcfR('{input}', verbose = F)

vcf.info <- as_tibble(snp.vcf@fix) %>%
  select(ID, CHROM, POS, REF, ALT, INFO) %>%
  separate(INFO,
           .$INFO[[1]] %>% str_remove_all('=[0-9.]+') %>% str_split(';') %>% unlist(),
           sep=';[^=]+=', convert = T) %>%
  mutate(#Chr = str_remove(CHROM, 'chr') %>% factor(levels=c(1:22, 'X', 'Y')),
         GC = str_remove(GC, 'GC='),
         POS = as.numeric(POS),
         alleles = paste0(REF, ',', ALT),
         A_freq = (2*N_AA + N_AB) / (2*(N_AA + N_AB + N_BB)),
         #This is B_allele frequency, or 'population frequency B-allele'. PennCNV needs this
         PFB = (2*N_BB + N_AB) / (2*(N_AA + N_AB + N_BB)),
         #This is AF/Alternate allele frequency. Some bcftool modules need this
         Alt_freq = ifelse(ALLELE_B == 1, PFB, A_freq)
         )
# ensure that CHROM uses NCBI/ENSEMBL styles
styles <- genomeStyles()
vcf.info$Chr <- factor(mapSeqlevels(vcf.info$CHROM %>% as.character(), 'NCBI'),
        levels = styles$Homo_sapiens$NCBI)

vcf.info %>%
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
    container:
        "docker://genomicslab/penncnv"
    shell:
        """
        if [[ -f /home/user/PennCNV/gc_file/{GENOME}.gc5Base.txt ]]; then
            ln -s /home/user/PennCNV/gc_file/{GENOME}.gc5Base.txt {params.download_path}/{GENOME}.gc5Base.txt
        else
            gunzip -c /home/user/PennCNV/gc_file/{GENOME}.gc5Base.txt.gz > {params.download_path}/{GENOME}.gc5Base.txt 
        fi
        /home/user/PennCNV/cal_gc_snp.pl {params.download_path}/{GENOME}.gc5Base.txt {input} -out {output}
        """


rule create_array_info_file:
    input:
        pfb=config["penncnv_pfb_file"],
        chromInfo=config["genomeInfo_file"],
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

# The 'getChromInfoFromUCSC' function breaks if UCSC does updates faster than Bioconductor
# -> use the chromInfo.txt file from UCSC instead
chromInfo <- read_tsv('{input.chromInfo}', show_col_types = FALSE) %>%
    mutate(chr = factor(chr, levels=unique(chr))) %>%
    group_by(chr) %>% summarise(seqlength = unique(size))
ginfo <- Seqinfo(seqnames=chromInfo$chr %>% as.character(),
                 seqlengths=chromInfo$seqlength,
                 isCircular= chromInfo$chr == 'chrMT',
                 genome='{params.genome}')

# Read pfb to get probe positions, pfb uses NCBI chromosome style
array <- read_tsv('{input.pfb}', show_col_types = FALSE) %>%
  as_granges(seqnames = Chr, start = Position, width = 1) %>%
  #Need to filter duplicate positions, as they screw up the density
  reduce_ranges(n_probes = n())
array <- renameSeqlevels(array,
                         mapSeqlevels(seqnames(array) %>% as.character(), 'UCSC') %>%
                           str_replace('chrM', 'chrMT') %>% setNames(seqnames(array) %>% as.character())
)
seqinfo(array) <- ginfo

# Get gaps
gaps <- complement_ranges(array) %>%
  mutate(chrom_start = start == 1,
         chrom_end = end == seqlengths(.)[as.character(seqnames)],
         ) %>%
  filter(!chrom_end & !chrom_start) %>%
  mutate(width_mean_sd = mean(width(.))+sd(width(.)),
         gap_min_size = ifelse('{params.min_gap_size}' == 'auto-array',
                               unique(width_mean_sd),
                               as.numeric('{params.min_gap_size}'))) %>%
  filter(width(.) >= gap_min_size)
write_bed(gaps, '{output.gaps}')

# calculate probe density 
# by tiling the whole genome we have even windows, and mean & median are similar
density <- chromInfo %>% as_granges(seqnames = chr, width = seqlength, start = 1) %>%
    tile_ranges({params.density_windows}) %>%
    mutate(score = count_overlaps(., array) / width(.))
write_bed(density, '{output.density}')

EOF
"""


rule gencode_v45_gtf_download:
    output:
        config["genome_gtf_file"],
    # Source gtf GRCh38:
    params:
        ftp_base="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/",
        release_path=(
            "gencode.v45.basic.annotation.gtf.gz"
            if GENOME == "hg38"
            else "GRCh37_mapping/gencode.v45lift37.basic.annotation.gtf.gz"
        ),
    shell:
        """
        wget {params.ftp_base}/{params.release_path} -O {output} 2> /dev/null
        # gunzip {output}.gz
        """


# #FIXME: the gencode fa.gz files are gzip not bgzip compressed
# rule gencode_v45_genomeFasta_download:
#     output: config['genome_fasta_file']
#     # Source fasta GRCh38:
#     params:
#         ftp_base = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/",
#         release_path = "GRCh38.p14.genome.fa.gz" if GENOME == "hg38" else "GRCh37_mapping/GRCh37.primary_assembly.genome.fa.gz"
#     shell:
#         """
#         wget {params.ftp_base}/{params.release_path} -O {output} 2> /dev/null
#         """


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
        config["genomeInfo_file"],
    conda:
        "../envs/general-R.yaml"
    shell:
        """
Rscript - << 'EOF'
suppressMessages(library(tidyverse))
suppressMessages(library(GenomeInfoDb))

styles <- genomeStyles()
styles_MT <- styles$Homo_sapiens$UCSC %>% str_replace('chrM', 'chrMT')

chrominfo <- read_tsv("{input.chrominfo}", col_names = c("chr", "size", "url"), show_col_types = FALSE) %>% 
    select(-url) %>% 
    filter(chr != 'chrMT') %>%
    mutate(chr = str_replace(chr, 'chrM$', 'chrMT')) %>%
    filter(chr %in% styles_MT)
# The 'acen' labelled cytobands are the centromere areas (but much larger than what the centromere file gives)
chrominfo <- left_join(chrominfo,
                read_tsv("{input.cytobands}", show_col_types = FALSE, 
                    col_names = c("chr", "band_start", "band_end", "band_name", "band_staining")),
                by = "chr"
                ) %>%
             mutate(
                centromer = band_staining == "acen",
                chr = factor(chr, styles_MT)
             ) %>%
             arrange(chr, band_start) %>%
             # The hg19 cytoband file does not include chrMT, so no band_start there
             filter(!is.na(band_start) | chr == 'chrMT')
write_tsv(chrominfo, "{output}")
EOF
        """


rule download_vep_fasta:
    output:
        os.path.join(
            config["vep_fasta_path"],
            "homo_sapiens",
            "112_{genome}",
            "Homo_sapiens.{genome}.dna.toplevel.fa.gz",
        ),
    conda:
        importlib.resources.files(STEM_CNV_CHECK).joinpath(
            "envs", "vep-annotation.yaml"
        )
    params:
        fasta_path=config["vep_fasta_path"],
    shell:
        "vep_install -a f -s homo_sapiens -y {wildcards.genome} -c {params.fasta_path}"


rule download_vep_cache:
    output:
        done=os.path.join(config["vep_cache_path"], ".{genome}.done"),
        folder=directory(
            os.path.join(config["vep_cache_path"], "homo_sapiens", "112_{genome}")
        ),
    conda:
        importlib.resources.files(STEM_CNV_CHECK).joinpath(
            "envs", "vep-annotation.yaml"
        )
    params:
        cache_path=config["vep_cache_path"],
    shell:
        # -r {params.cache_path}/plugins; check if needed or if it follows -c
        # -a f > automatically get ensembl genome fasta; can replace fasta & gtf
        "vep_install -a cp -g DosageSensitivity -s homo_sapiens -y {wildcards.genome} -c {params.cache_path} --CONVERT && touch {output.done}"
