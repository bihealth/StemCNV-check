# Redirect warnings & errors to snakemake logging, save R environment if debugging
source(file.path(snakemake@config$snakedir, 'scripts/common.R'))

library(tidyverse)
library(DNAcopy)
library(plyranges)
library(vcfR)
library(GenomeInfoDb)

snakemake@source('R/helper_functions.R')
snakemake@source('R/vcf_io_functions.R')
snakemake@source('R/CNV_preprocess_functions.R')


CBS_LRR_segmentation <- function(tb, CBS_config, sex, target_style, sample_id = 'test') {
   
    # CBS code
    cna.basic <- CNA(tb$LRR, tb$seqnames, tb$start, data.type = 'logratio', sampleid = sample_id)
    cna.basic.smoothed <- smooth.CNA(cna.basic) #Note: can we also apply GC-wave correction here?
    cna.basic.smoothed.segmented <- segment(
        cna.basic.smoothed, 
        # "minimum number of markers for a changed segment"
        # -> does not actually remove segments <X probes, capped at 2-5
        min.width = 5, 
        undo.splits = 'sdundo', 
        undo.SD = CBS_config$undo.SD.val
    )

    # Need sex chroms in matching style to set proper CNs
    sex_chroms <- get_sex_chroms(target_style)
    
    # Re-formatting
    tb <- segments.summary(cna.basic.smoothed.segmented) %>%
  		dplyr::rename(seqnames = chrom, start = loc.start, end = loc.end, n_probes = num.mark) %>%
        fix_CHROM_format(target_style) %>%
        mutate(
            sample_id = sample_id,
            width = end - start + 1,
            CN = case_when(
                #Male is default CN=1 on X & Y, also unqiue cutoffs
                sex == 'm' & seqnames %in% sex_chroms & seg.median < CBS_config$LRR.male.XorY.loss       ~ 0,
                sex == 'm' & seqnames %in% sex_chroms & seg.median > CBS_config$LRR.male.XorY.gain.large ~ 3,
                sex == 'm' & seqnames %in% sex_chroms & seg.median > CBS_config$LRR.male.XorY.gain       ~ 2,
                sex == 'm' & seqnames %in% sex_chroms                                                    ~ 1,
                # Unique cutoffs for female X (behaves different from autosome)
                sex == 'f' & seqnames == sex_chroms[1] & seg.median < CBS_config$LRR.female.XX.loss       ~ 0,
                sex == 'f' & seqnames == sex_chroms[1] & seg.median < CBS_config$LRR.female.X.loss        ~ 1,
                sex == 'f' & seqnames == sex_chroms[1] & seg.median > CBS_config$LRR.female.X.gain.large  ~ 4,
                sex == 'f' & seqnames == sex_chroms[1] & seg.median > CBS_config$LRR.female.X.gain        ~ 3,
                sex == 'f' & seqnames == sex_chroms[1]                                                    ~ 2,
                # Default cutoffs
                seg.median < CBS_config$LRR.loss.large                                                    ~ 0,
                seg.median < CBS_config$LRR.loss                                                          ~ 1,
                seg.median > CBS_config$LRR.gain.large                                                    ~ 4,
                seg.median > CBS_config$LRR.gain                                                          ~ 3,
                TRUE ~ 2,
                .default = 2
            ),
            CNV_type = case_when(# Male is default CN=1 on X & Y
                #FIXME (future): technically DUP should ONLY be used for CN=3 (or 2 on male XY)
                sex == 'm' & seqnames %in% sex_chroms & CN == 2 ~ 'DUP',
                sex == 'm' & seqnames %in% sex_chroms & CN == 1 ~ NA,
                CN > 2                                          ~ 'DUP',
                CN < 2                                          ~ 'DEL',
                .default = NA
            ),
            CNV_caller = 'CBS',
            ID = paste(CNV_caller, CNV_type, seqnames, start, end, sep='_'),
        ) %>%
        filter(!is.na(CNV_type) & !is.na(seqnames)) %>%
        select(seqnames, start, end, width, ID, CNV_caller, CNV_type, CN, sample_id)

    if (sex == 'f') {
        tb <- filter(tb, seqnames != sex_chroms[2])
    }
    
    return(tb)
    
}


get_CBS_CNV_vcf <- function(input_vcf, out_vcf, config, sample_id = 'test') {
    # Get settings
    sample_sex <- get_sample_info(sample_id, "sex", config)
    tool_config <- config$settings$CBS
    defined_labels <- get_defined_labels(config)
    
    # Get SNP data
    snp_vcf <- read.vcfR(input_vcf, verbose = F) 
    snp_vcf_meta <- snp_vcf@meta
    snp_vcf_gr <- parse_snp_vcf(snp_vcf)
    target_style <- get_target_chrom_style(config, snp_vcf_gr)
    
    # Run CBS
    cnv_gr <- CBS_LRR_segmentation(as_tibble(snp_vcf_gr), tool_config, sample_sex, target_style, sample_id) %>%
        as_granges()
    
    #make sure that snp_vcf & cnv_vcf use the same (& intended) chrom style
    cnv_gr <- fix_CHROM_format(cnv_gr, target_style)
    snp_vcf_gr <- fix_CHROM_format(snp_vcf_gr, target_style)
    
    # preprocess (merge, filter, SNP counts)
    cnvs <- apply_preprocessing(cnv_gr, snp_vcf_gr, tool_config) %>%
        get_median_LRR(snp_vcf_gr) %>%
        as_tibble()
    
    # Write VCF
    filtersettings <- tool_config$probe_filter_settings
    if (filtersettings == '_default_') {
        filtersettings <- config$settings$default_probe_filter_set
    } 
    vcf_info_text <- paste(
        '##CBS=R-DNAcopy LRR segmentation: CNA/smooth.CNA/segment',
        str_glue('undo.SD={tool_config$undo.SD.val} min.width=5'),
        str_glue('CN_LRR_thresholds_autosomes: CN0={tool_config$LRR.loss.large} CN1={tool_config$LRR.loss}'),
        str_glue('CN3={tool_config$LRR.gain} CN4={tool_config$LRR.gain.large}'),
        str_glue('CN_LRR_thresholds_female_X: CN0={tool_config$LRR.female.XX.loss} CN1={tool_config$LRR.female.X.loss}'),
        str_glue('CN3={tool_config$LRR.female.X.gain} CN4={tool_config$LRR.female.X.fain.large}'), 
        str_glue('CN_LRR_thresholds_male_XY: CN0={tool_config$LRR.male.XorY.loss}'),
        str_glue('CN2={tool_config$LRR.male.XorY.gain} CN3={tool_config$LRR.male.XorY.gain.large}'),
        str_glue('StemCNV-check_array_probe_filtering="{filtersettings}"')
    )
    write_cnv_vcf(
        cnvs,
        out_vcf,
        sample_sex,
        'CBS',
        config,
        defined_labels,
        snp_vcf_meta,
        vcf_info_text,
        target_style
    )
}


get_CBS_CNV_vcf(
    input_vcf = snakemake@input[['vcf']],
    out_vcf = snakemake@output[['vcf']],
    config = snakemake@config,
    sample_id = snakemake@wildcards[['sample_id']]
)