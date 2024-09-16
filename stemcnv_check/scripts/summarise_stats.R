# Redirect warnings & errors to snakemake logging, save R environment if debugging
source(file.path(snakemake@config$snakedir, 'scripts/common.R'))

library(tidyverse)
library(furrr)
library(readxl)
library(writexl)

snakemake@source('R/helper_functions.R')
snakemake@source('R/vcf_io_functions.R')
snakemake@source('R/processCNVs_annotate_check-score.R')

plan(multisession, workers = snakemake@threads)

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
    s %in% config$evaluation_settings$summary_stat_warning_levels$last_level_critical,
    'critical',
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

get_summary_overview_table <- function(gencall_stats, sample_SNP_gr, SNP_distance_to_reference, sample_CNV_data, config) {
    
    sample_id <- unique(gencall_stats$sample_id)

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
            sample_id = sample_id,
            SNP_distance_to_reference = SNP_distance_to_reference
        ),
        get_call_stats(sample_CNV_data)
    )
    
    callrate_warnings <- config$evaluation_settings$summary_stat_warning_levels$call_rate
    sample_sex <- get_sample_info(sample_id, 'sex', config$sample_table)
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
            SNPs_post_filter = NA_character_
        ),
        full_join(qc_measure_list[[3]],qc_measure_list[[4]]) %>%
            mutate(across(2:(length(qc_measure_list[[4]])+1), ~apply_greq_th(., cur_column(), config)))
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

get_call_stats <- function(gr.or.tb, name_addition = NA) {
    
    tb <- gr.or.tb %>%
        as_tibble() %>%
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

sample_GT_distances <- function(sample_SNP_gr, extra_snp_files, use_filter=TRUE) {
    # Reading the SNP vcf files is the slowest thing here
    # > parallelized via furrr:future_map to speed it up, esp. for large number of files
    # FIXME (future): maybe don't parse the actual VCF here and use a faster table reader (readr / data.table) instead? 
    future_map(
        extra_snp_files,
        parse_snp_vcf,
        format_fields = c('GT'),
        apply_filter = use_filter,
        .options = furrr_options(seed = TRUE)
    ) %>%
        bind_ranges( sample_SNP_gr) %>%
        select(ID, sample_id, GT) %>%
        as_tibble() %>%
        mutate(GT = ifelse(GT == './.', NA, str_count(GT, '1'))) %>%
        pivot_wider(names_from = sample_id, values_from = GT, values_fill = NA) %>%
        filter(if_all(everything(), ~!is.na(.))) %>%
        dplyr::select(-seqnames, -start, -end, -strand, -width, -ID) %>%
        t() %>%
        dist(method = 'manhattan') %>%
        as.matrix() %>% 
        as.data.frame() %>%
        as_tibble(rownames = 'sample_distance_to')
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
    
    sample_id <- unique(gencall_stats$sample_id)
    ref_id <- get_sample_info(sample_id, 'ref_id', config$sample_table)
    SNP_distance_to_reference <- ifelse(
        is.na(ref_id),
        NA_real_,
        SNP_GT_distances %>%
            filter(sample_distance_to == ref_id) %>%
            pull(!!sample_id)
    )
    
    summary_table_sample <- get_summary_overview_table(
        gencall_stats,
        sample_SNP_gr,
        SNP_distance_to_reference,
        sample_CNV_gr,
        config
    )
    if (!is.null(summary_excel_ref)) {
        summary_table_sample <- summary_table_sample %>%
            full_join(
                read_excel(summary_excel_ref, sheet = 'summary_stats') %>%
                    rename_with(~ str_replace(., 'sample', 'reference')) %>%
                    mutate(across(matches('critical|reportable'), ~ NA_character_)),
                by = 'Description'
            )
    }
        
    tool_stats <- unsplit_merged_CNV_callers(sample_CNV_gr) %>%
        split(.$CNV_caller) %>%
        imap(function (gr, name) {
            tb1 <- gr %>%
                annotate_call.label(config$evaluation_settings$CNV_call_categorisation) %>%
                get_call_stats()
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
        list(gencall_stats = gencall_stats %>% tr_sample_tb()),
        set_names(tool_stats, paste0(names(tool_stats), '_stats')),
        list(PennCNV_QC_info = parse_penncnv_logs(penncnv_log_files)),
        list(SNP_GT_distances = SNP_GT_distances)
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
    write_xlsx(snakemake@output[['xlsx']])