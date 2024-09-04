# Redirect all output to snakemake logging
log <- file(snakemake@log[['err']], 'wt')
sink(log, append = T)
sink(log, append = T, type = 'message')

library(tidyverse)
library(readxl)
library(writexl)

snakemake@source('R/helper_functions.R')
snakemake@source('R/vcf_io_functions.R')

tr_sample_tb <- function(tb) {
  tr <- t(tb )
  colnames(tr) <- c('Value')
  as_tibble(tr, rownames = 'Description')
}

#TODO: this stay in the report 
# get_sample_stats <- function(sample_id, config) {
#     sampletable <- read_sampletable(config$sample_table)
#     use_ids <- c(sample_id, get_sample_info(sample_id, 'ref_id', sampletable)) %>% na.omit()
#     
#     out <- lapply(use_ids, function(id) {
#         sampletable %>%
#             filter(Sample_ID == id) %>%
#             dplyr::select(
#                 Sample_ID, Sex, Reference_Sample, config$evaluation_settings$summary_stats$sample.info.extra.cols
#             ) %>%
#             mutate(`Analysis run date` = Sys.Date()) %>%
#             tr_sample_tb() %>%
#             dplyr::rename(sample_value = Value)
#     }) %>%
#         full_join() %>%
# }     

get_gencall_stats <- function(gencall_stat_file) {
    read_tsv(gencall_stat_file, show_col_types = FALSE, 
             col_types = cols(computed_gender = col_character(), sentrix_barcode = col_character())
    ) %>%
        dplyr::select(-(5:7)) %>%
        mutate(gtc = str_remove(gtc, '.gencall.gtc')) %>%
        dplyr::rename(sample_id = gtc)
}


get_summary_overview_table <- function(gencall_stats, sample_SNP_gr, SNP_GT_distances, sample_CNV_data, config) {
    
    not_used_filter <- c(
        ifelse(config$evaluation_settings$CNV_call_categorisation$impact.score.critical == 'NA', 'critical', 'dummy'),
        ifelse(config$evaluation_settings$CNV_call_categorisation$impact.score.reportable == 'NA', 'reportable', 'dummy')
    ) %>% unique %>% paste(collapse = '|')
    
    
    if ('Call_Label' %!in% colnames(mcols(sample_CNV_data))) {
        sample_CNV_data <- add_call_labels(sample_CNV_data, config)
    }

    qc_measure_list <- list(
        gencall_stats %>%
            dplyr::select(sample_id, call_rate, computed_gender) %>%
            mutate(call_rate = round(call_rate, 3)),
        sample_SNP_gr %>%
            as_tibble() %>%
            group_by(sample_id) %>%
            summarise(
                SNPs_post_filter = (100 * sum(FILTER == 'PASS', na.rm=T) / unique(gencall_stats$number_snps) ) %>%
                    format(digits = 2, nsmall=2) %>% paste('%')
            ),
        tibble(
            sample_id = gencall_stats$sample_id,
            SNP_distance_to_reference = SNP_GT_distances[1]
        ),
        get_call_stats(sample_CNV_data, filter_regex = not_used_filter)
    )
    
    # Each measure used for QC is categorised as 'normal', 'abnormal', or 'problematic'
    # CURRENT: 'OK', 'unusual', 'warning', 'critical'
    # Additonally the last level can eb set to 'critical' instead of 'problematic'  
    # The report color codes 'normal' as green, 'abnormal' as yellow, 'problematic' as orange, and 'critical' as red
    get_last_level <- function(s) ifelse(
        s %in% config$evaluation_settings$summary_stat_warning_levels$last_level_critical,
        'critical',
        'warning'
    )

    apply_greq_th <- function(datacol, measure){
        #FIXME (future): Error or default for what happens when not defined?
        warning_levels <- config$evaluation_settings$summary_stat_warning_levels[[measure]]
        case_when(
            is.na(datacol)                 ~ NA_character_,
            datacol >= warning_levels[[2]] ~ get_last_level(measure),
            datacol >= warning_levels[[1]] ~ 'warning',
            TRUE                           ~ 'OK',
        )
    }
    
    callrate_warnings <- config$evaluation_settings$summary_stat_warning_levels$call_rate
    sample_sex <- get_sample_info(sample_id, 'sex', config$sampletable)
    qc_eval_list <- list(
        mutate(
            qc_measure_list[[1]],
            call_rate = case_when(
                call_rate < callrate_warnings[[2]]  ~ get_last_level('call_rate'),
                call_rate < callrate_warnings[[1]]  ~ 'warning',
                TRUE ~ 'OK'
            ),
            computed_gender = ifelse(tolower(computed_gender) != sample_sex, get_last_level('computed_gender'), 'OK')
        ),
        mutate(
            qc_measure_list[[2]],
            SNPs_post_filter = NA_character_
        ),
        full_join(qc_measure_list[[3]],qc_measure_list[[4]]) %>%
            mutate(across(2:(length(qc_measure_list[[4]])+1), ~apply_greq_th(., cur_column())))
    )

    
    full_join(        
        purrr::reduce(qc_measure_list, full_join, by = 'sample_id') %>%
            tr_sample_tb() %>%
            dplyr::rename(sample_value = Value), 
        purrr::reduce(qc_eval_list, full_join, by = 'sample_id') %>%
            tr_sample_tb() %>%
            dplyr::rename(sample_eval = Value)
    )

}


#TODO any reason to NOT do this in the preprocessing ? (beside the report not getting it that way yet)
add_call_labels <- function(gr, config) {
    
    impact.score.critical <- config$evaluation_settings$CNV_call_categorisation$impact.score.critical
    critical_excl_regex <- config$evaluation_settings$CNV_call_categorisation$filters.exclude.critical %>%
        c('dummy') %>% paste(collapse = '|')
    impact.score.reportable <- config$evaluation_settings$CNV_call_categorisation$impact.score.reportable
    reportable_excl_regex <- config$evaluation_settings$CNV_call_categorisation$filters.exclude.reportable %>%
        c('dummy') %>% paste(collapse = '|')
    
    gr %>%
        as_tibble() %>%
        mutate(
            Call_Label = case_when(
                !is.na(reference_coverage)     				                                        ~ 'Reference genotype',
                Check_Score >= impact.score.critical & !str_detect(FILTER, critical_excl_regex)     ~ 'Critical',
                Check_Score >= impact.score.reportable & !str_detect(FILTER, reportable_excl_regex) ~ 'Reportable',
                TRUE                         				                                        ~ '...'
            )
        )
}


unsplit_merged_CNV_callers <- function(cnv_gr) {
    
    unmerged_calls <- cnv_gr %>%
        filter(CNV_caller != 'StemCNV-check')
    
    merged_calls <- cnv_gr %>%
        filter(CNV_caller == 'StemCNV-check') %>%
        as_tibble() %>%
        separate_rows(initial_call_details, sep = '\\|') %>%
        mutate(
            CNV_caller = str_extract(initial_call_details, '^[^_]+'),
            start = str_extract(initial_call_details, '(?<=_)[0-9]+(?=-)') %>% as.integer(),
            end = str_extract(initial_call_details, '(?<=-)[0-9]+(?=_CN)') %>% as.integer(),
            CN = str_extract(initial_call_details, '(?<=CN)[0-9]+(?=_cov)') %>% as.integer(),
            # overlap_merged_call = str_extract(initial_call_details, '(?<=cov)[0-9]+\\.[0-9]+')
        ) %>%
        select(-width) %>%
        as_granges()
    
    bind_ranges(merged_calls, unmerged_calls)

}


get_call_stats <- function(gr.or.tb, name_addition = NA, filter_regex = NA) {
    
    tb <- gr.or.tb %>%
        group_by(sample_id) %>%
        mutate(
            loss_gain_log2ratio = log2(sum(CNV_type == 'gain') / sum(CNV_type == 'loss')) %>% round(digits = 2),
            loss_gain_log2ratio = ifelse(is.infinite(loss_gain_log2ratio) | is.nan(loss_gain_log2ratio), NA, loss_gain_log2ratio),
            highlight_calls = sum(!is.na(Highlight)),
            high_impact_calls = sum(!is.na(HighImpact)),
            CNV_type = ifelse(CNV_type == 'LOH', 'LOH', 'CNV'),
        ) %>%
        group_by(sample_id, CNV_type, loss_gain_log2ratio) %>%
        summarise(
            total_calls = dplyr::n(),
            reportable_calls = sum(Call_Label == 'Reportable'),
            critical_calls = sum(Call_Label == 'Critical')
        ) %>%
        pivot_wider(names_from = CNV_type, values_from = matches('(total|reportable|critical)_calls')) 

    if (!is.na(filter_regex)) {
        tb <- tb %>%
            dplyr::select(-matches(filter_regex))
    }
    
    if (!is.na(name_addition)) {
        tb <- tb %>%
            rename_with(~paste(., name_addition, sep = '_'), !matches('sample_id'))
    }
    
    tb
}


# Extract qc summary stats from PennCNV log files
parse_penncnv_logs <- function(penncnv_logfiles) {
    
    lapply(penncnv_logfiles, function(fname) {
        chrs <- basename(fname) %>% str_extract('^[^\\.]+')
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
        dplyr::rename(` ` = Name)
}


sample_GT_distances <- function(sample_SNP_gr, extra_snp_files, use_filter=TRUE) {
    
    extra_snp_files %>% # as.list() %>%
        # FIXME (future): This, sadly, is pretty slow
        # maybe parallelize via furrr:future_map instead?
        # -> set threads to n_callers in snakemake (or default 2?)
        # -> plan(multisession, workers = snakemake@threads)
        map(
            parse_snp_vcf,
            format_fields = c('GT'),
            apply_filter = use_filter
        ) %>%
        bind_ranges( sample_SNP_gr) %>%
        select(ID, sample_id, GT) %>%
        as_tibble() %>%
        mutate(GT = ifelse(GT == './.', NA, str_count(GT, '1'))) %>%
        pivot_wider(names_from = sample_id, values_from = GT, values_fill = NA) %>%
        filter(if_all(everything(), ~!is.na(.))) %>%
        dplyr::select(-seqnames, -start, -end, -strand, -width, -ID) %>%
        t() %>%
        dist(method = 'manhattan')
}

collect_summary_stats <- function(
    gencall_stat_file,
    SNP_vcf_file,
    #TODO: SNP_analysis_file,
    CNV_vcf_file,
    penncnv_log_files,
    extra_snp_files,
    config,
    summary_excel_ref = NULL
) {
    
    sample_SNP_gr <- parse_snp_vcf(SNP_vcf_file, format_fields = c('GT'), apply_filter = FALSE)
    sample_CNV_gr <- parse_cnv_vcf(CNV_vcf_file, apply_filter = FALSE) 
    
    
    gencall_stats <- get_gencall_stats(gencall_stat_file)
    
    SNP_GT_distances <- sample_GT_distances(
        sample_SNP_gr, 
        extra_snp_files, 
        use_filter = config$evaluation_settings$SNP_clustering$`filter-settings` != 'none'
    )
    # SNP_GT_distances %>% as.matrix() %>% as.data.frame() %>% as.dist()
    
    
    summary_table_sample <- get_summary_overview_table(
        gencall_stats,
        sample_SNP_gr,
        SNP_GT_distances,
        sample_CNV_gr,
        config
    )
    if (!is.null(summary_excel_ref)) {
        summary_table_sample <- summary_table_sample %>%
            full_join(
                read_excel(summary_excel_ref, sheet = 'sample_info') %>%
                    rename_with(~ str_replace(., 'sample', 'reference')) %>%
                    mutate(across(matches('critical|reportable'), ~ NA_character_)),
                by = 'Description'
            )
    }
    
    not_used_filter <- c(
        ifelse(config$evaluation_settings$CNV_call_categorisation$impact.score.critical == 'NA', 'critical', 'dummy'),
        ifelse(config$evaluation_settings$CNV_call_categorisation$impact.score.reportable == 'NA', 'reportable', 'dummy')
    ) %>% unique %>% paste(collapse = '|')
    
    tool_stats <- unsplit_merged_CNV_callers(sample_CNV_gr) %>%
        split(.$CNV_caller) %>%
        imap(function (gr, name) {
            tb1 <- gr %>%
                add_call_labels(config) %>%
                get_call_stats(
                    # name_addition = name,
                    filter_regex = not_used_filter
                )
            tb2 <- tb1 %>%
                 mutate(across(2:(ncol(tb2)-1), ~apply_greq_th(., cur_column())))
            
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
                            rename_with(~ str_replace(., 'sample', 'reference')) %>%
                            mutate(across(matches('critical|reportable'), ~ NA_character_)),
                        by = 'Description'
                    )
            }
            out
                
        })
    
    # write to multiple excel sheets
    # - summary stat table
    # - tool stat tables
    # - PennCNV log stats (too tricky to merge with tool stat table)
    #FIXME (future): move this to somewhere else?
    # - SNP_GT_distances 
    
    return( c(
        list(summary_stats = summary_table_sample),
        set_names(tool_stats, paste0(names(tool_stats), '_stats')),
        list(PennCNV_QC_info = parse_penncnv_logs(penncnv_log_files)),
        list(SNP_GT_distances = as.matrix(SNP_GT_distances) %>% as.data.frame() %>% as_tibble(rownames = 'sample_distance_to'))
    ))
               
    
}


collect_summary_stats(
    gencall_stat_file = snakemake@input[['gencall_stats']],
    SNP_vcf_file = snakemake@input[['snp_vcf']],
    #TODO: SNP_analysis_file,
    CNV_vcf_file = snakemake@input[['cnv_vcf']],
    penncnv_log_files = snakemake@input[['penncnv_logs']],
    extra_snp_files = snakemake@input[['extra_snp_files']],
    config = snakemake@config,
    summary_excel_ref = snakemake@input[['summary_excel_ref']]
) %>% 
    write_xlsx(snakemake@output)