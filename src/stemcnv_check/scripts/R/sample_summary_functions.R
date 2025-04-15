suppressMessages(library(tidyverse))
suppressMessages(library(yaml))

# Helper function, could be defined elsewhere?
get_defined_labels <- function(config) {
    label_def_file <- file.path(config$snakedir, 'control_files', 'label_name_definitions.yaml')
    return(yaml.load_file(label_def_file))
}

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
        mutate(gtc = str_remove(gtc, '.gencall.[a-zA-Z0-9]+.gtc')) %>%
        dplyr::rename(sample_id = gtc)
}

# sample levels include names and colors (4 each) for QC measures and are defined in the 'label_name_definitions.yaml' file
# The last level is only used for certain measures (based on config)
# Default names should be: 'OK', 'unusual', 'warning' / 'high concern'
# Default colors are: green, yellow, orange / red
apply_measure_th <- function(
    datacol, measure, sample_levels, config, comp = 'greater_equal', th_override = NULL
){
    stopifnot(comp %in% c('greater_equal', 'smaller', 'not_equal'))
    comp_f <- list(
        'greater_equal' = function(x,y) {x >= y},
        'smaller'       = function(x,y) {x < y},
        'not_equal'     = function(x,y) {x != y}
    ) [[comp]]
    
    #FIXME (future): Error or default for what happens when not defined?
    measure_thresholds <- config$evaluation_settings$summary_stat_warning_levels[[measure]]
    last_level <-  measure %in% config$evaluation_settings$summary_stat_warning_levels$use_last_level

    if (!is.null(th_override)) {
        measure_thresholds <- th_override
    }    
    
    case_when(
        is.na(datacol)                                        ~ NA_character_,
        comp_f(datacol, measure_thresholds[[2]]) & last_level ~ sample_levels[[4]],
        comp_f(datacol, measure_thresholds[[2]])              ~ sample_levels[[3]],
        comp_f(datacol, measure_thresholds[[1]])              ~ sample_levels[[2]],
        TRUE                                                  ~ sample_levels[[1]],
    )
}


get_summary_overview_table <- function(gencall_stats, snp_qc_details, sample_CNV_gr, config, defined_labels) {
    
    sample_id <- unique(gencall_stats$sample_id)

    qc_measure_list <- list(
        gencall_stats %>%
            dplyr::select(sample_id, call_rate, computed_gender) %>%
            mutate(call_rate = round(call_rate, 3)),
        snp_qc_details %>%
            dplyr::select(
                sample_id, SNPs_post_filter, SNP_pairwise_distance_to_reference, reportable_SNVs, critical_SNVs
            ) %>%
            unique(),
        get_call_stats(sample_CNV_gr, config$evaluation_settings$summary_stat_warning_levels$call_count_excl_labels)
    )
    
    sample_levels <- names(get_defined_labels(config)$sample_labels)    
    sample_sex <- get_sample_info(sample_id, 'sex', config)
    
    qc_eval_list <- list(
        mutate(
            qc_measure_list[[1]],
            call_rate = apply_measure_th(
                call_rate, 'call_rate', sample_levels, config, 'smaller'
            ),
            computed_gender = apply_measure_th(
                tolower(computed_gender), 'computed_gender', sample_levels, config,
                comp = 'not_equal', th_override = list(sample_sex, sample_sex)
            )
        ),
        mutate(
            qc_measure_list[[2]],
            SNPs_post_filter = NA_character_,
            SNP_pairwise_distance_to_reference = apply_measure_th(
                SNP_pairwise_distance_to_reference, 'SNP_pairwise_distance_to_reference', sample_levels, config
            ),
            reportable_SNVs = ifelse(
                is.na(SNP_pairwise_distance_to_reference),
                NA_character_,
                apply_measure_th(reportable_SNVs, 'reportable_SNVs', sample_levels, config)
            ),
            critical_SNVs = ifelse(
                is.na(SNP_pairwise_distance_to_reference),
                NA_character_,
                apply_measure_th(critical_SNVs, 'critical_SNVs', sample_levels, config)
            )
        ),
        qc_measure_list[[3]] %>%
            mutate(across(where(is.numeric), ~apply_measure_th(., cur_column(), sample_levels, config)))
    )
    
    tb <- full_join(        
        purrr::reduce(qc_measure_list, full_join, by = 'sample_id') %>%
            tr_sample_tb() %>%
            dplyr::rename(sample_value = Value), 
        purrr::reduce(qc_eval_list, full_join, by = 'sample_id') %>%
            tr_sample_tb() %>%
            dplyr::rename(sample_eval = Value)
    )
    # Sort the rows based on the order of the defined labels
    tb %>%
        mutate(Description = factor(Description, levels = c('sample_id', defined_labels$sample_qc_measures))) %>%
        arrange(Description)
}


get_call_stats <- function(gr.or.tb, call_count_excl_labels = list(), name_addition = NA) {
    
    stopifnot(
        is(gr.or.tb, 'GRanges') | is(gr.or.tb, 'GRangesList') | is(gr.or.tb, 'data.frame'),
        is(call_count_excl_labels, 'list') | is(call_count_excl_labels, 'character')
    )
    
    tb <- gr.or.tb %>%
        as_tibble() %>%
        filter(Call_label %!in% call_count_excl_labels) %>%
        group_by(sample_id) %>%
        mutate(
            loss_gain_log2ratio = log2(sum(CNV_type == 'gain') / sum(CNV_type == 'loss')) %>% round(digits = 2),
            loss_gain_log2ratio = ifelse(is.infinite(loss_gain_log2ratio) | is.nan(loss_gain_log2ratio), NA, loss_gain_log2ratio),
            cancer_gene_calls = sum(!is.na(cancer_gene)),
            stemcell_hotspot_calls = sum(!is.na(stemcell_hotspot)),
            CNV_type = ifelse(CNV_type == 'LOH', 'LOH', 'CNV'),
        ) %>%
        group_by(sample_id, CNV_type, loss_gain_log2ratio) %>%
        # Note: the call labels are now config definable!
        # >> Names should NOT be hard-coded here?!
        # >> either define via config or defined_labels
        summarise(
            total_calls = dplyr::n(),
            reportable_calls = sum(Call_label == 'Reportable de-novo', na.rm = TRUE),
            critical_calls = sum(Call_label == 'Critical de-novo', na.rm = TRUE)
        ) %>%
        pivot_wider(names_from = CNV_type, values_from = matches('(total|reportable|critical)_calls'))
    
    if (!is.na(name_addition)) {
        tb <- tb %>%
            rename_with(~paste(., name_addition, sep = '_'), !matches('sample_id'))
    }
    
    tb
}
