# Redirect warnings & errors to snakemake logging, save R environment if debugging
source(file.path(snakemake@config$snakedir, 'scripts/common.R'))

library(tidyverse)
library(readxl)
library(writexl)

snakemake@source('R/helper_functions.R')
snakemake@source('R/vcf_io_functions.R')
snakemake@source('R/CNV_annotation_functions.R')
snakemake@source('R/sample_summary_functions.R')


collect_summary_stats <- function(
    sample_id,
    gencall_stat_file,
    SNP_analysis_file,
    CNV_vcf_file,
    penncnv_log_files,
    config,
    summary_excel_ref = NULL
) {
    
    defined_labels <- get_defined_labels(config)
    
    sample_CNV_gr <- parse_cnv_vcf(CNV_vcf_file, apply_filter = FALSE)
    
    gencall_stats <- get_gencall_stats(gencall_stat_file)
    
    stopifnot(sample_id == unique(gencall_stats$sample_id))
    ref_id <- get_sample_info(sample_id, 'ref_id', config)
    
    snp_qc_details <- read_excel(SNP_analysis_file, sheet = 'SNP_QC_details')
    
    summary_table_sample <- get_summary_overview_table(
        gencall_stats,
        snp_qc_details,
        sample_CNV_gr,
        config,
        defined_labels
    )
    if (!is.na(ref_id)) {
        summary_table_sample <- summary_table_sample %>%
            full_join(
                read_excel(summary_excel_ref, sheet = 'summary_stats') %>%
                    select(-contains('reference')) %>%
                    rename_with(~ str_replace(., 'sample', 'reference')) %>%
                    filter(!str_detect(Description, 'critical|reportable')),
                by = 'Description'
            )
    }
    
    min.ref.coverage <- config$settings$CNV_processing$call_processing$min.reciprocal.coverage.with.ref
    sample_levels <- names(defined_labels$sample_labels)
    tool_stats <- split_merged_CNV_callers(sample_CNV_gr, defined_labels) %>%
        split(.$CNV_caller) %>%
        imap(function (gr, name) {
            tb1 <- gr %>%
                mutate(reference_overlap = ifelse(
                    is.na(reference_coverage), FALSE, reference_coverage >= min.ref.coverage
                )) %>%
                annotate_call.label(config$evaluation_settings$CNV_call_labels) %>%
                get_call_stats(config$evaluation_settings$summary_stat_warning_levels$call_count_excl_labels)
            tb2 <- tb1 %>%
                 mutate(across(2:(ncol(tb1)-1), ~apply_measure_th(., cur_column(), sample_levels, config)))
            
            out <- full_join(
                tb1 %>%
                    tr_sample_tb() %>%
                    dplyr::rename(sample_value = Value), 
                tb2 %>%
                    tr_sample_tb() %>%
                    dplyr::rename(sample_eval = Value)
            )
                
            if (!is.null(summary_excel_ref)) {
                out <- out %>%
                    full_join(
                        read_excel(summary_excel_ref, sheet = paste0(name, '_stats')) %>%
                            select(-contains('reference')) %>%
                            rename_with(~ str_replace(., 'sample', 'reference')) %>%
                            filter(!str_detect(Description, 'critical|reportable')),
                        by = 'Description'
                    )
            }
            out
                
        })
    
    # write to multiple excel sheets
    # - summary stat table
    # - tool stat tables
    # - PennCNV log stats (too tricky to merge with tool stat table)
    
    return( c(
        list(summary_stats = summary_table_sample),
        list(gencall_stats = gencall_stats %>% tr_sample_tb()),
        set_names(tool_stats, paste0(names(tool_stats), '_stats')),
        list(PennCNV_QC_info = parse_penncnv_logs(penncnv_log_files))
    ))
               
    
}


collect_summary_stats(
    sample_id = snakemake@wildcards[['sample_id']],
    gencall_stat_file = snakemake@input[['gencall_stats']],
    SNP_analysis_file = snakemake@input[['snv_analysis']],
    CNV_vcf_file = snakemake@input[['cnv_vcf']],
    penncnv_log_files = snakemake@input[['penncnv_logs']],
    config = snakemake@config,
    summary_excel_ref = snakemake@input[['summary_excel_ref']]
) %>% 
    write_xlsx(snakemake@output[['xlsx']])