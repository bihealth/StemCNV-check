# Redirect warnings & errors to snakemake logging, save R environment if debugging
source(file.path(snakemake@config$snakedir, 'scripts/common.R'))

library(tidyverse)
#library(furrr)
library(readxl)
library(writexl)

snakemake@source('R/helper_functions.R')
snakemake@source('R/vcf_io_functions.R')
snakemake@source('R/CNV_annotation_functions.R')

# plan(multisession, workers = snakemake@threads)

tr_sample_tb <- function(tb) {
  tr <- t(tb )
  colnames(tr) <- c('Value')
  as_tibble(tr, rownames = 'Description')
}

get_gencall_stats <- function(gencall_stat_file) {
    read_tsv(gencall_stat_file, show_col_types = FALSE, 
             col_types = cols(computed_gender = col_character(), sentrix_barcode = col_character())
    ) %>%
        dplyr::select(-(5:7)) %>%
        mutate(gtc = str_remove(gtc, '.gencall.gtc')) %>%
        dplyr::rename(sample_id = gtc)
}

# Each measure used for QC is categorised as 'OK', 'unusual', or 'warning'
# Additonally the last level can eb set to 'critical' instead of 'warning'  
# The report color codes 'OK' as green, 'unusual' as yellow, 'warning' as orange, and 'critical' as red
get_last_level <- function(s, config) ifelse(
    s %in% config$evaluation_settings$summary_stat_warning_levels$last_level_red,
    'high concern',
    'warning'
)

apply_greq_th <- function(datacol, measure, config){
    #FIXME (future): Error or default for what happens when not defined?
    warning_levels <- config$evaluation_settings$summary_stat_warning_levels[[measure]]
    case_when(
        is.na(datacol)                 ~ NA_character_,
        datacol >= warning_levels[[2]] ~ get_last_level(measure, config),
        datacol >= warning_levels[[1]] ~ 'unusual',
        TRUE                           ~ 'OK',
    )
}

get_summary_overview_table <- function(gencall_stats, snp_qc_details, sample_CNV_gr, config) {
    
    sample_id <- unique(gencall_stats$sample_id)

    qc_measure_list <- list(
        gencall_stats %>%
            dplyr::select(sample_id, call_rate, computed_gender) %>%
            mutate(call_rate = round(call_rate, 3)),
        snp_qc_details %>%
            dplyr::select(sample_id, SNPs_post_filter, SNP_pairwise_distance_to_reference, critical_snvs) %>%
            unique(),
        get_call_stats(sample_CNV_gr, config$evaluation_settings$CNV_call_categorisation$call_count_excl_filters)
    )
    
    callrate_warnings <- config$evaluation_settings$summary_stat_warning_levels$call_rate
    sample_sex <- get_sample_info(sample_id, 'sex', config)
    qc_eval_list <- list(
        mutate(
            qc_measure_list[[1]],
            call_rate = case_when(
                call_rate < callrate_warnings[[2]]  ~ get_last_level('call_rate', config),
                call_rate < callrate_warnings[[1]]  ~ 'unusual',
                TRUE ~ 'OK'
            ),
            computed_gender = ifelse(tolower(computed_gender) != sample_sex, get_last_level('computed_gender', config), 'OK')
        ),
        mutate(
            qc_measure_list[[2]],
            SNPs_post_filter = NA_character_,
            SNP_pairwise_distance_to_reference = apply_greq_th(SNP_pairwise_distance_to_reference, 'SNP_pairwise_distance_to_reference', config),
            critical_snvs = ifelse(
                is.na(SNP_pairwise_distance_to_reference),
                NA_character_,
                apply_greq_th(critical_snvs, 'critical_snvs', config)
            )
        ),
        # full_join(qc_measure_list[[3]],)
        qc_measure_list[[3]] %>%
            mutate(across(where(is.numeric), ~apply_greq_th(., cur_column(), config)))
    )
    
    tb <- full_join(        
        purrr::reduce(qc_measure_list, full_join, by = 'sample_id') %>%
            tr_sample_tb() %>%
            dplyr::rename(sample_value = Value), 
        purrr::reduce(qc_eval_list, full_join, by = 'sample_id') %>%
            tr_sample_tb() %>%
            dplyr::rename(sample_eval = Value)
    )
    # Sort the "critical" columns together
    bind_rows(
        tb %>% filter(Description != 'critical_snvs'),
        tb %>% filter(Description == 'critical_snvs')
    )
}

get_call_stats <- function(gr.or.tb, call_count_excl_filters = list(), name_addition = NA) {
    
    call_filter_regex <- ifelse(
        is.null(call_count_excl_filters) || length(call_count_excl_filters) == 0,
        'dummy',
        call_count_excl_filters %>% paste(collapse = '|')
    )
    
    tb <- gr.or.tb %>%
        as_tibble() %>%
        filter(is.na(FILTER) | !str_detect(FILTER, call_filter_regex)) %>%
        group_by(sample_id) %>%
        mutate(
            loss_gain_log2ratio = log2(sum(CNV_type == 'gain') / sum(CNV_type == 'loss')) %>% round(digits = 2),
            loss_gain_log2ratio = ifelse(is.infinite(loss_gain_log2ratio) | is.nan(loss_gain_log2ratio), NA, loss_gain_log2ratio),
            cancer_gene_calls = sum(!is.na(cancer_gene)),
            stemcell_hotspot_calls = sum(!is.na(stemcell_hotspot)),
            CNV_type = ifelse(CNV_type == 'LOH', 'LOH', 'CNV'),
        ) %>%
        group_by(sample_id, CNV_type, loss_gain_log2ratio) %>%
        summarise(
            total_calls = dplyr::n(),
            reportable_calls = sum(Call_label == 'Reportable', na.rm = TRUE),
            critical_calls = sum(Call_label == 'Critical', na.rm = TRUE)
        ) %>%
        pivot_wider(names_from = CNV_type, values_from = matches('(total|reportable|critical)_calls'))
    
    if (!is.na(name_addition)) {
        tb <- tb %>%
            rename_with(~paste(., name_addition, sep = '_'), !matches('sample_id'))
    }
    
    tb
}


# Extract qc summary stats from PennCNV log files
parse_penncnv_logs <- function(penncnv_log_files) {
    
    lapply(penncnv_log_files, function(fname) {
        chrs <- basename(fname) %>% str_extract('(?<=PennCNV.)[^.]+')
        lines <- readLines(fname)
        # Extract from summary line
        tb <- lines %>%
          str_subset('NOTICE: quality summary') %>%
                str_remove('.*: ') %>%
                str_split(' ') %>%
                unlist() %>%
                as_tibble_col(column_name = 'dummy') %>%
                separate(dummy, c('Name', 'value'), sep = '=') %>%
                filter(Name %!in% c('WF', 'GCWF')) %>%
                mutate(Name = str_remove(Name, '[XY]'))
        # Get wave correction
        tb <- bind_rows(tb, lines %>%
                str_subset('Adjusting LRR by GC model') %>%
                str_remove('.*: ') %>%
                str_split(', ') %>%
                unlist() %>%
                as_tibble_col(column_name = 'dummy') %>%
                separate(dummy, c('Name', 'value'), sep = ' changes from ') %>%
                mutate(Name = paste(Name, '(adjusted)'))
        )
        # Get median correction BAF & LRR (auto only)
        if (chrs == 'auto') {
            tb[str_detect(tb$Name, 'median' ),] <- lines %>%
                    str_subset('Median-adjusting') %>%
                    str_replace('.*(LRR|BAF).*( -?[0-9]\\.[0-9]{4})', '\\1_median;0-adjusted by\\2') %>%
                    as_tibble() %>%
                    separate(value, c('Name', 'value'), sep = ';')
        }
    
        mutate(tb, chr = ifelse(chrs == 'auto', 'chr1:22', chrs))
    }) %>%
        bind_rows() %>%
        pivot_wider(names_from = chr, values_from = value) %>%
        dplyr::rename(Description = Name)
}


collect_summary_stats <- function(
    sample_id,
    gencall_stat_file,
    SNP_analysis_file,
    CNV_vcf_file,
    penncnv_log_files,
    config,
    summary_excel_ref = NULL
) {
    
    sample_CNV_gr <- parse_cnv_vcf(CNV_vcf_file, apply_filter = FALSE)
    
    gencall_stats <- get_gencall_stats(gencall_stat_file)
    
    stopifnot(sample_id == unique(gencall_stats$sample_id))
    ref_id <- get_sample_info(sample_id, 'ref_id', config)
    
    snp_qc_details <- read_excel(SNP_analysis_file, sheet = 'SNP_QC_details')
    
    summary_table_sample <- get_summary_overview_table(
        gencall_stats,
        snp_qc_details,
        sample_CNV_gr,
        config
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
        
    tool_stats <- unsplit_merged_CNV_callers(sample_CNV_gr) %>%
        split(.$CNV_caller) %>%
        imap(function (gr, name) {
            tb1 <- gr %>%
                annotate_call.label(config$evaluation_settings$CNV_call_categorisation) %>%
                get_call_stats(config$evaluation_settings$summary_stat_warning_levels$call_count_excl_filters)
            tb2 <- tb1 %>%
                 mutate(across(2:(ncol(tb1)-1), ~apply_greq_th(., cur_column(), config)))
            
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