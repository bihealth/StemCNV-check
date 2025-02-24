# Redirect warnings & errors to snakemake logging, save R environment if debugging
source(file.path(snakemake@config$snakedir, 'scripts/common.R'))

library(tidyverse)
library(plyranges)
library(writexl)

snakemake@source('R/helper_functions.R')
snakemake@source('R/vcf_io_functions.R')

cnv_collate_call_selection <- snakemake@params$cnv_collate_call_selection

tb <- lapply(snakemake@input$cnv_calls, parse_cnv_vcf) %>%
    bind_ranges() %>%
    as_tibble() %>%
    dplyr::select(-strand) %>%
    dplyr::rename(chromosome = seqnames, vcf_filter = FILTER)

# Change default column order
tb <- tb %>% 
    # sample_id, coords, ID, FILTER, CNV_type, Call_label, Check_Score, Prec. Estimate, Hotspots, 
    select(
        sample_id, 1:6, CNV_type, 12, 10, 11, 13:16,
        any_of(colnames(tb))
    )

if (!is.null(cnv_collate_call_selection$whitelist_call_label)) {
    tb <- tb %>% filter(Call_label %in% output_filters$whitelist_call_label)
}

if (!is.null(cnv_collate_call_selection$blacklist_call_label)) {
    tb <- tb %>% filter(Call_label %!in% output_filters$blacklist_call_label)
}

if (snakemake@params$output_format == 'xlsx') {
    write_xlsx(tb, snakemake@output[[1]])
} else {
    # tsv
    write_tsv(tb, snakemake@output[[1]])
}