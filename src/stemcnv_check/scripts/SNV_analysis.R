# Redirect warnings & errors to snakemake logging, save R environment if debugging
source(file.path(snakemake@config$snakedir, 'scripts/common.R'))

library(tidyverse)
library(writexl)
library(furrr)
library(plyranges)
library(vcfR)
library(GenomeInfoDb)

snakemake@source('R/helper_functions.R')
snakemake@source('R/vcf_io_functions.R')
snakemake@source('R/hotspot_functions.R')
snakemake@source('R/CNV_annotation_functions.R')
snakemake@source('R/snv_analysis_functions.R')

#multicore support
plan(multisession, workers = snakemake@threads)


run_snp_analysis <- function(
    sample_id,
    sample_SNP_vcf_file,
    ref_SNP_vcf_file,
    extra_snp_files,
    config,
    gtf_file,
    ginfo_file
) {
    # Internal tables
    defined_labels <- get_defined_labels(config)
    SNV_hotspot_table <- load_hotspot_table(config, 'snv_hotspot')
    # sample info
    ref_id <- get_sample_info(sample_id, 'ref_id', config) 
    
    ## Load SNP data for sample
    sample_SNP_gr <- parse_snp_vcf(
        sample_SNP_vcf_file,
        info_fields = c('GenTrain_Score', 'ANN'),
        format_fields = c('GT', 'IGC'),
        apply_filter = FALSE
    )
    target_chrom_style <- get_target_chrom_style(config, sample_SNP_gr)
    sample_SNP_gr <- fix_CHROM_format(sample_SNP_gr, target_chrom_style)
    # and for ref (if it exists)
    if(!is.na(ref_id)) {
        ref_SNP_gr <- parse_snp_vcf(
            ref_SNP_vcf_file,
            format_fields = c('GT', 'IGC'),
            apply_filter = FALSE
        ) %>%
            fix_CHROM_format(target_chrom_style)
    } else {
        ref_SNP_gr <- GRanges(
            REF = character(),
            ALT = character(),
            IGC = integer(),
            GT = character()
        )
    }
    # Prepare internal table of SNVs: 
    # - annotate SNVs with sample ROI
    # - extract mehari annotations into separate columns
    sample_SNV_tb <- get_sample_SNV_tb(sample_SNP_gr, sample_id, SNV_hotspot_table, gtf_file, ginfo_file, target_chrom_style, config)
    # Prepare table for xlsx output
    # - compare with reference GT
    # - check which critical reasons apply (ROI, SNV-hotspot, HIGH impact, ..)
    # - label SNVs (critical, non-critical, etc.)
    SNV_table <- get_SNV_table(sample_SNV_tb, ref_SNP_gr, SNV_hotspot_table, config, defined_labels)
    
    # Calculate sample distance matrix
    SNP_GT_distances <- sample_GT_distances(
        sample_SNP_gr, 
        ref_SNP_gr,
        extra_snp_files,
        ref_SNP_vcf_file,
        target_chrom_style,
        use_filter = config$settings$SNV_analysis$`probe_filter_settings` != 'none'
    )
    
    # Collect output tables for xlsx file
    return(list(
        'SNV_table' = SNV_table,
        'SNV_hotspot_coverage' = get_SNV_hotspot_coverage(sample_SNV_tb, SNV_hotspot_table),
        'SNP_GT_distances' = SNP_GT_distances,
        'SNP_QC_details' = get_SNV_QC_table(
            sample_id, sample_SNV_tb, ref_SNP_gr, SNV_table, 
            use_filter = config$settings$SNV_analysis$`probe_filter_settings` != 'none'
        )
    ))
    
}

run_snp_analysis(
    sample_id = snakemake@wildcards[['sample_id']],
    sample_SNP_vcf_file = snakemake@input[['snp_vcf']],
    ref_SNP_vcf_file = snakemake@input[['ref_snp_vcf']],
    extra_snp_files = snakemake@input[['extra_snp_files']],
    config = snakemake@config,
    gtf_file = snakemake@params[['gtf_file']],
    ginfo_file = snakemake@params[['ginfo_file']]
) %>% 
    write_xlsx(snakemake@output[['xlsx']])

