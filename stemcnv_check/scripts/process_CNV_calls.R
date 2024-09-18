# Redirect warnings & errors to snakemake logging, save R environment if debugging
source(file.path(snakemake@config$snakedir, 'scripts/common.R'))

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

snp_vcf <- read.vcfR(snakemake@input$snp_vcf, verbose = F) 
snp_vcf_meta <- snp_vcf@meta
snp_vcf_gr <- parse_snp_vcf(snp_vcf)
target_chrom_style <- get_target_chrom_style(config, snp_vcf_gr)

snp_vcf_gr <- snp_vcf_gr %>%
    fix_CHROM_format(target_chrom_style)
gr <- load_cnv_callers(snakemake@input$cnv_calls) %>%
    # should not be necessary here, but we do it anyway 
    fix_CHROM_format(target_chrom_style)

# Tool overlapping
combined_tools_sample <- combine_CNV_callers(gr, processing_config, snp_vcf_gr) 
## }

## <- function() {
if (length(snakemake@input$ref_data) > 0) {
    
	combined_tools_ref <- parse_cnv_vcf(snakemake@input$ref_data) %>%
		# These cols will interefere
		select(-reference_coverage, -Genes)
    
	
	cnvs <- annotate_reference_overlap(combined_tools_sample, combined_tools_ref,
                                       processing_config$min.reciprocal.coverage.with.ref)
} else {
	cnvs <- combined_tools_sample %>%
        mutate(
            reference_overlap = FALSE,
            reference_coverage = NA_real_,
            reference_caller = NA_character_
        )
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

gr_genes <- load_gtf_data(config, target_chrom_style)
gr_info  <- load_genomeInfo(config, target_chrom_style)

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
    plyranges::select(-LRR) %>%
    get_median_LRR(snp_vcf_gr) %>%
	annotate_impact_lists(high_impact_gr, 'high_impact') %>%
	annotate_impact_lists(highlight_gr, 'highlight') %>%
	annotate_roi(sample_id, sampletable, gr_genes, gr_info) %>%
	annotate_gaps(config$static_data$array_gaps,
                  processing_config$min.perc.gap_area, 
                  processing_config$gap_area.uniq_probes.rel) %>%
	annotate_high_density(config$static_data$array_density, processing_config$density.quantile.cutoff) %>%
	annotate_gene_overlaps(gr_genes) %>%
    as_tibble() %>%
	annotate_cnv.check.score(high_impact_gr, highlight_gr, check_scores) %>%
	annotate_precision.estimates(size_categories, precision_estimates) %>%
    annotate_call.label(config$evaluation_settings$CNV_call_categorisation)


# Also directly write out a cnv vcf
combined_calls_to_vcf <- function(cnv_tb, vcf_out, sample_sex, processing_config, vcf_meta, target_style) {
    
    tb <- cnv_tb %>%
        as_tibble() %>%
        rowwise() %>%
        mutate(
            CNV_type = ifelse(CNV_type == 'LOH', 'CNV:LOH', CNV_type),
            CNV_caller = ifelse(length(CNV_caller) > 1, 'StemCNV-check', unlist(CNV_caller)),
            FILTER = ifelse(is.na(FILTER), 'PASS', FILTER),
        )
    
    filtersettings <- processing_config$`filter-settings`
    if (filtersettings == '_default_') {
        filtersettings <- config$settings$`default-filter-settings`
    }    

    header <- c(
        fix_header_lines(vcf_meta, 'fileformat|contig|BPM=|EGT=|CSV=', target_style),
        static_cnv_vcf_header(processing_config, extra_annotation = TRUE),
        '##ALT=<ID=CNV:LOH,Description="Loss of heterozygosity, same as run of homozygosity">',
        #FIXME (future): maybe also keep the PennCNV & CBS lines? (doesn't seem to be 100% standard though)
        str_glue(
            '##StemCNV-check process_CNV_calls',
            'tool.overlap.greatest.call.min.perc={processing_config$tool.overlap.greatest.call.min.perc}',
            'tool.overlap.min.cov.sum.perc={processing_config$tool.overlap.min.cov.sum.perc}'
        )
    )
    
    cnv_vcf <- new(
        "vcfR",
        meta = header,
        fix = get_fix_section(tb),
        gt = get_gt_section(tb, sample_sex)
    )
    
    write.vcf(cnv_vcf, vcf_out)
    
}

cnvs %>%
    combined_calls_to_vcf(
        snakemake@output$vcf,
        get_sample_info(sample_id, 'sex', sampletable),
        processing_config,
        snp_vcf_meta,
        target_chrom_style
    )
