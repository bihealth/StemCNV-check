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
snakemake@source('R/CNV_preprocess_functions.R')
snakemake@source('R/CNV_comparison_functions.R')
snakemake@source('R/CNV_annotation_functions.R')
snakemake@source('R/hotspot_functions.R')

##################
# Variable setup #
##################

config <- snakemake@config
sample_id <- snakemake@wildcards$sample_id
sampletable <- read_sampletable(config$sample_table, config$column_remove_regex)
sample_sex <- get_sample_info(sample_id, 'sex', config, sampletable)

# Load data & set CHROM Style
processing_config <- config$settings$CNV_processing$call_processing
defined_labels <- get_defined_labels(config)

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
combined_tools_sample <- combine_CNV_callers(gr, processing_config, snp_vcf_gr, defined_labels) 
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

stemcell_hotspots_tb <- load_hotspot_table(config, 'stemcell_hotspot')
cancer_genes_tb <- load_hotspot_table(config, 'cancer_gene')
dosage_sensitive_gene_tb <- get_dosage_sensivity_tb(
    snakemake@params$dosage_file,
    config
)
roi_tb <- get_roi_tb(sample_id, sampletable, config)
hotspot_genes <- bind_rows(
    stemcell_hotspots_tb,
    cancer_genes_tb,
    dosage_sensitive_gene_tb,
    roi_tb
) %>%
    filter(mapping == 'gene_name') %>%
    pull(hotspot) %>%
    unique()

gr_genes <- load_gtf_data(snakemake@params$gtf_file, config, target_chrom_style, include_hotspot_genes = hotspot_genes)
gr_info  <- load_genomeInfo(snakemake@params$ginfo_file, config, target_chrom_style)

stemcell_hotspots_gr <- parse_hotspot_table(stemcell_hotspots_tb, gr_genes, gr_info)
cancer_genes_gr <- parse_hotspot_table(cancer_genes_tb, gr_genes, gr_info)
dosage_sensitive_gene_gr <- parse_hotspot_table(dosage_sensitive_gene_tb, gr_genes, gr_info)

array <- sampletable %>%
    filter(Sample_ID == sample_id) %>%
    pull('Array_Name')
gap_file <- config$array_definition[[array]]$array_gaps_file
density_file <- config$array_definition[[array]]$array_density_file

cnvs <- cnvs %>%
    plyranges::select(-LRR) %>%
    get_median_LRR(snp_vcf_gr) %>%
	annotate_impact_lists(stemcell_hotspots_gr, 'stemcell_hotspot') %>%
    annotate_impact_lists(dosage_sensitive_gene_gr, 'dosage_sensitive_gene') %>%
	annotate_impact_lists(cancer_genes_gr, 'cancer_gene') %>%
	annotate_roi(roi_tb, gr_genes, gr_info, config) %>%
	annotate_gaps(
        gap_file,
        processing_config$min.perc.gap_area, 
        processing_config$gap_area.uniq_probes.rel,
        target_chrom_style
    ) %>%
	annotate_high_density(
        density_file,
        processing_config$density.quantile.cutoff,
        target_chrom_style
    ) %>%
	annotate_gene_overlaps(gr_genes) %>%
    as_tibble() %>%
	annotate_cnv.check.score(stemcell_hotspots_gr, dosage_sensitive_gene_gr, cancer_genes_gr, check_scores, sample_sex) %>%
	annotate_precision.estimates(size_categories, precision_estimates) %>%
    annotate_call.label(config$evaluation_settings$CNV_call_labels)


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
    
    filtersettings <- processing_config$`probe_filter_settings`
    if (filtersettings == '_default_') {
        filtersettings <- config$settings$default_probe_filter_set
    }    

    header <- c(
        fix_header_lines(vcf_meta, 'fileformat|contig|BPM=|EGT=|CSV=', target_style),
        static_cnv_vcf_header(processing_config, extra_annotation = TRUE, fullconfig = config),
        '##ALT=<ID=CNV:LOH,Description="Loss of heterozygosity, same as run of homozygosity">',
        #FIXME (future): maybe also keep the PennCNV & CBS lines? (doesn't seem to be 100% standard though)
        str_glue(
            '##StemCNV-check process_CNV_calls',
            'tool.overlap.greatest.call.min.perc={processing_config$tool.overlap.greatest.call.min.perc}',
            'tool.overlap.min.cov.sum.perc={processing_config$tool.overlap.min.cov.sum.perc}'
        )
    )
    
    fix <- get_fix_section(tb)
    gt <- get_gt_section(tb, sample_id, sample_sex, target_style)
    # write.vcf does not work on empty vcfR objects
    if (nrow(tb) == 0) {
        vcf_file <- out_vcf %>% str_replace('.gz$', '')
        cat(header, file = vcf_file, sep = '\n')
        cat(
            paste0(
                '#', paste(c(colnames(fix), colnames(gt)), collapse = '\t')
            ),
            file = vcf_file, sep = '\n', append = T
        )
        R.utils::gzip(vcf_file)
    } else {
        cnv_vcf <- new(
            "vcfR",
            meta = header,
            fix = fix,
            gt = gt
        )
        write.vcf(cnv_vcf, out_vcf)
    }
}

cnvs %>%
    combined_calls_to_vcf(
        snakemake@output$vcf,
        get_sample_info(sample_id, 'sex', snakemake$config, sampletable),
        processing_config,
        snp_vcf_meta,
        target_chrom_style
    )
