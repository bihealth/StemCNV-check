# Redirect all output to snakemake logging
log <- file(snakemake@log[['err']], 'wt')
sink(log, append = T)
sink(log, append = T, type = 'message')

library(plyranges)
library(GenomicRanges)
library(tidyverse)
library(vcfR)
library(yaml)
`%!in%` <- Negate(`%in%`)
options(dplyr.summarise.inform = FALSE)

snakemake@source('R/helper_functions.R')
snakemake@source('R/vcf_io_functions.R')
snakemake@source('R/preprocess_CNV_functions.R')

snakemake@source('R/processCNVs_combine_CNV_callers.R')
snakemake@source('R/processCNVs_annotate_reference_overlap.R')
snakemake@source('R/processCNVs_annotate_impact_lists.R')
snakemake@source('R/processCNVs_annotate_array_features.R')
snakemake@source('R/processCNVs_annotate_check-score.R')

##################
# Variable setup #
##################

config <- snakemake@config
sample_id <- snakemake@wildcards$sample_id
sampletable <- read_sampletable(config$sample_table)

# Load data & set CHROM Style
processing_config <- config$settings$CNV_processing$call_processing


## <- function() {
target_chrom_style <- config$settings$vcf_output$chrom_style
if (target_chrom_style == '__keep__') { target_style <- seqlevelsStyle(snp_vcf_gr) %>% head(1) 
} else { target_style <- target_chrom_style }

snp_vcf <- parse_snp_vcf(snakemake@input$snp_vcf) %>%
    fix_CHROM_format(target_chrom_style)
gr <- load_cnv_callers(snakemake@input$cnv_calls) %>%
    # should not be necessary here, but we do it anyway 
    fix_CHROM_format(target_chrom_style)

# Tool overlapping
combined_tools_sample <- combine_CNV_callers(gr, processing_config, snp_vcf)
## }

## <- function() {
if (length(snakemake@input$ref_data) > 0) {
    
	combined_tools_ref <- load_preprocessed_cnvs(snakemake@input$ref_data) %>%
		# These cols will interefere
		dplyr::select(-width, -reference_overlap, -reference_coverage, -reference_caller,
									-n_genes, -overlapping_genes) %>%
		filter(caller_merging_state != 'pre-overlap') %>%
		as_granges()
	
	cnvs <- annotate_reference_overlap(combined_tools_sample, combined_tools_ref,
                                       processing_config$min.reciprocal.coverage.with.ref)
} else {
	cnvs <- combined_tools_sample
}
## }


# # Gaps
# gap_area.uniq_probes.rel <- config$settings$CNV_processing$call_processing$gap_area.uniq_probes.rel
# min.perc.gap_area <- config$settings$CNV_processing$call_processing$min.perc.gap_area
# # High Density
# density.quantile.cutoff <- config$settings$CNV_processing$call_processing$density.quantile.cutoff



#Scoring
check_scores <- config$settings$CNV_processing$Check_score_values
#Precision esitmation
size_categories <- config$settings$CNV_processing$Precision$size_categories
precision_estimates<- config$settings$CNV_processing$Precision$estimate_values

gr_genes <- load_gtf_data(config)
gr_info  <- load_genomeInfo(config)

HI_file <- config$settings$CNV_processing$gene_overlap$high_impact_list %>%
	str_replace('__inbuilt__', config$snakedir)
if (!is.null(HI_file)) {
	high_impact_gr <- read_tsv(HI_file, show_col_types = FALSE) %>%
		parse_hotspot_table(gr_genes, gr_info)
} else {
	high_impact_gr <- GRanges()
}
HL_file <- config$settings$CNV_processing$gene_overlap$highlight_list %>%
	str_replace('__inbuilt__', config$snakedir)
if (!is.null(HL_file)) {
	highlight_gr <- read_tsv(HL_file, show_col_types = FALSE) %>%
		parse_hotspot_table(gr_genes, gr_info)
} else {
	highlight_gr <- GRanges()
}



cnvs <- cnvs %>%
	annotate_impact_lists(high_impact_gr, 'high_impact') %>%
	annotate_impact_lists(highlight_gr, 'highlight') %>%
	annotate_roi(sample_id, sampletable, gr_genes, gr_info) %>%
	annotate_gaps(config$static_data$array_gaps,
                  processing_config$min.perc.gap_area, 
                  processing_config$gap_area.uniq_probes.rel) %>%
	annotate_high_density(config$static_data$array_density, processing_config$density.quantile.cutoff) %>%
	finalise_gr_to_tb(gr_genes) %>%
	annotate_cnv.check.score(high_impact_gr, highlight_gr, check_scores) %>%
	annotate_precision.estimates(size_categories, precision_estimates)

#TODO: write cnv vcf

cnvs.tb <- cnvs %>%
	rowwise() %>%
	mutate(across(one_of(get_list_cols()), ~paste(., collapse=';')))
write_tsv(cnvs.tb, snakemake@output$tsv)
