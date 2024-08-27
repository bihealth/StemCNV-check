# Redirect all output to snakemake logging
log <- file(snakemake@log[['err']], 'wt')
sink(log, append = T)
sink(log, append = T, type = 'message')

library(tidyverse)
library(DNAcopy)
library(plyranges)
library(vcfR)
library(GenomeInfoDb)

snakemake@source('R/helper_functions.R')
snakemake@source('R/vcf_io_functions.R')
snakemake@source('R/preprocess_CNV_functions.R')


CBS_LRR_segmentation <- function(tb, CBS_config, sex, sample_id = 'test') {
   
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
    sex_chroms <- get_sex_chroms(tb)
    
    # Re-formatting
    tb <- segments.summary(cna.basic.smoothed.segmented) %>%
  		dplyr::rename(seqnames = chrom, start = loc.start, end = loc.end, n_probes = num.mark) %>%
        mutate(
            seqnames = str_remove(seqnames, '^(chr|ch)'),
            sample_id = sample_id,
            width = end - start + 1,
            # Chr = paste0('chr', Chr),
            # Chr = factor(Chr, levels = c(paste0('chr', 1:22), 'chrX', 'chrY')),
            # snp.density = n_snp_probes / length * 1e6,
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
                #TODO: technically DUP should ONLY be used for CN=3 (or 2 on male XY)
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
    sample_sex <- get_sample_info(sample_id, "sex", config$sample_table)
    tool_config <- config$settings$CBS
    # Get SNP data
    snp_vcf <- read.vcfR(input_vcf, verbose = F) 
    snp_vcf_meta <- snp_vcf@meta
    snp_vcf_gr <- parse_snp_vcf(snp_vcf)
    # Run CBS
    cnv_gr <- CBS_LRR_segmentation(as_tibble(snp_vcf_gr), tool_config, sample_sex, sample_id) %>%
        as_granges()
    
    #make sure that snp_vcf & cnv_vcf use the same (& intended) chrom style
    target_style <- get_target_chrom_style(config, snp_vcf_gr)
    cnv_gr <- fix_CHROM_format(cnv_gr, target_style)
    snp_vcf_gr <- fix_CHROM_format(snp_vcf_gr, target_style)
    
    # preprocess (merge, filter, SNP counts)
    cnvs <- apply_preprocessing(cnv_gr, snp_vcf_gr, tool_config) %>%
        get_median_LRR(snp_vcf_gr) %>%
        as_tibble()
    # Generate VCF
    filtersettings <- tool_config$`filter-settings`
    if (filtersettings == '__default__') {
        filtersettings <- config$settings$`default-filter-settings`
    }
    header <- c(
        fix_header_lines(snp_vcf_meta, 'fileformat|contig|BPM=|EGT=|CSV=', target_style),
        static_cnv_vcf_header(tool_config),
        paste(
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
    )
    
    cnv_vcf <- new(
        "vcfR",
        meta = header,
        fix = get_fix_section(cnvs),
        gt = get_gt_section(cnvs, sample_sex)
    )
    
    write.vcf(cnv_vcf, out_vcf)
    
}


get_CBS_CNV_vcf(
    snakemake@input[['vcf']],
    snakemake@output[['vcf']],
    snakemake@config,
    snakemake@wildcards[['sample_id']]
)