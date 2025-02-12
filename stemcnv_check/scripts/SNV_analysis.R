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
    
    sampletable <- read_sampletable(config$sample_table)
    ref_id <- get_sample_info(sample_id, 'ref_id', config) 
    
    ## Load SNP data from sample (& ref)
    sample_SNP_gr <- parse_snp_vcf(
        sample_SNP_vcf_file,
        info_fields = c('GenTrain_Score', 'ANN'),
        format_fields = c('GT', 'IGC'),
        apply_filter = FALSE
    ) 
    
    target_chrom_style <- get_target_chrom_style(config, sample_SNP_gr)
    sample_SNP_gr <- fix_CHROM_format(sample_SNP_gr, target_chrom_style)

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
    
    SNV_hotspot_table <- load_hotspot_table(config, 'snv_hotspot')
    
    roi_tb <- get_roi_tb(sample_id, sampletable, config)
    hotspot_genes <- c(
        roi_tb %>% filter(mapping == 'gene_name') %>% pull(hotspot),
        SNV_hotspot_table$gene_name
    ) %>% unique()
    gr_genes <- load_gtf_data(gtf_file, config, target_chrom_style, include_hotspot_genes = hotspot_genes)
    gr_info  <- load_genomeInfo(ginfo_file, config, target_chrom_style)
    
    # Sort, filter & Label annotated SNVs
    
    #Note: this should ideally be parsed from vcf header
    #Note: can we get rsIDs from here also? (they are not consistently used as SNP probe IDs)
    meahri_ann_header <- c(
        'Allele', 'Annotation', 'Annotation_Impact', 'gene_name', 'gene_id', 'Feature_Type', "Feature_ID",
        'Transcript_BioType', 'Rank', 'HGVS.c', 'HGVS.p', 'cDNA.pos', 'CDS.pos', 'AA.pos',        
        # cDNA.pos / cDNA.length | CDS.pos / CDS.length | AA.pos / AA.length |
        'Distance', 'Strand', 'errors' # ERRORS / WARNINGS / INFO
    )
    sample_SNV_tb <- sample_SNP_gr %>%
        annotate_roi(roi_tb, gr_genes, gr_info, config) %>%
        as_tibble() %>%
        separate(ANN, sep = '\\|', into = meahri_ann_header) %>%
        dplyr::rename(
            GenCall_Score = IGC,
            Transcript_ID = Feature_ID
        ) %>%
        mutate(
            GT = ifelse(is.na(GT), './.', GT)
        )
    
    # use filter here, yes or no?
    SNV_table <- get_SNV_table(sample_SNV_tb, ref_SNP_gr, SNV_hotspot_table, config)
    
    # Calculate sample distance matrix
    SNP_GT_distances <- sample_GT_distances(
        sample_SNP_gr, 
        ref_SNP_gr,
        extra_snp_files,
        ref_SNP_vcf_file,
        target_chrom_style,
        use_filter = config$settings$SNV_analysis$`filter-settings` != 'none'
    )
    
    # Collect output tables for xlsx file
    return(list(
        'SNV_table' = SNV_table,
        'SNV_hotspot_coverage' = get_SNV_hotspot_coverage(sample_SNV_tb, SNV_hotspot_table),
        'SNP_GT_distances' = SNP_GT_distances,
        'SNP_QC_details' = get_SNV_QC_table(
            sample_id, sample_SNV_tb, ref_SNP_gr, SNV_table, 
            use_filter = config$settings$SNV_analysis$`filter-settings` != 'none'
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

