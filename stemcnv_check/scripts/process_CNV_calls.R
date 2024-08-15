#! /usr/bin/Rscript
suppressMessages(library(argparse))

parser <- ArgumentParser(description="Collect data from CNV-calling tools and merge them into a single output, with added overlap stats")

parser$add_argument('data_path', type = 'character', help='Path to pipeline output data')
parser$add_argument('sample_id', type = 'character', help='Sample_ID to use')
parser$add_argument('configfile', type = 'character', help='Path to config file')
parser$add_argument('sampletablefile', type = 'character', help='Path to sampletable file')

parser$add_argument('-p', '--penncnv', action="store_true", help="Use PennCNV calls")
parser$add_argument('-c', '--cbs', action="store_true", help="Use CBS calls")

args <- parser$parse_args()

suppressMessages(library(plyranges))
suppressMessages(library(GenomicRanges))
suppressMessages(library(tidyverse))
suppressMessages(library(yaml))
`%!in%` <- Negate(`%in%`)
options(dplyr.summarise.inform = FALSE)

config <- read_yaml(args$configfile)

get_script_dir <- function() {
	match <- commandArgs(trailingOnly = FALSE) %>%
		str_subset('--file=') %>%
		str_remove('--file=')
	check <- sys.frames()[[1]]$ofile
	if (length(match) > 0) {
        # Rscript
		return(normalizePath(match) %>% dirname())
	} else if (!is.null(check)) {
		# source'd via R console
        return(normalizePath(check) %>% dirname())
	} else {
		# likely testing in IDE
		return ('stemcnv_check/scripts')
	}
}
source(file.path(get_script_dir(), 'R/R_io_functions.R'))
source(file.path(get_script_dir(), 'R/preprocess_CNV_functions.R'))
source(file.path(get_script_dir(), 'R/processCNVs_combine_CNV_callers.R'))
source(file.path(get_script_dir(), 'R/processCNVs_annotate_reference_overlap.R'))
source(file.path(get_script_dir(), 'R/processCNVs_annotate_impact_lists.R'))
source(file.path(get_script_dir(), 'R/processCNVs_annotate_array_features.R'))
source(file.path(get_script_dir(), 'R/processCNVs_annotate_check-score.R'))

##################
# Variable setup #
##################

datapath <- args$data_path
sample_id <- args$sample_id
sampletable <- read_tsv(args$sampletablefile, col_types = 'cccccc', comment = '#')

sex <- get_sample_info(sample_id, 'sex', sampletable)
ref_id <- get_sample_info(sample_id, 'ref_id', sampletable)
sex.ref <- get_sample_info(sample_id, 'sex.ref', sampletable)

# CNV callers
tool.order <- config$settings$CNV.calling.tools
# Call merging & pre-filtering
min.snp				<- config$settings$CNV_processing$call_processing$min.snp
min.length			<- config$settings$CNV_processing$call_processing$min.length
min.snp.density 	<- config$settings$CNV_processing$call_processing$min.snp.density
merge.distance  	<- config$settings$CNV_processing$call_processing$merge.distance
merge.before.filter <- config$settings$CNV_processing$call_processing$merge.before.filter
n_snp.filterset     <- ifelse(config$settings$CNV_processing$call_processing$`filter-settings` == '__default__',
                               config$settings$`default-filter-set`,
                               config$settings$CNV_processing$call_processing$`filter-settings`)
# Tool overlapping
tool.overlap.greatest.call.min.perc <- config$settings$CNV_processing$call_processing$tool.overlap.greatest.call.min.perc
tool.overlap.median.cov.perc <- config$settings$CNV_processing$call_processing$tool.overlap.median.cov.perc
# Compare to reference
min.reciprocal.coverage.with.ref <- config$settings$CNV_processing$call_processing$min.reciprocal.coverage.with.ref
# Gaps
gap_area.uniq_probes.rel <- config$settings$CNV_processing$call_processing$gap_area.uniq_probes.rel
min.perc.gap_area <- config$settings$CNV_processing$call_processing$min.perc.gap_area
# High Density
density.quantile.cutoff <- config$settings$CNV_processing$call_processing$density.quantile.cutoff


#Scoring
check_scores <- config$settings$CNV_processing$Check_score_values
#Precision esitmation
size_categories <- config$settings$CNV_processing$Precision$size_categories
precision_estimates<- config$settings$CNV_processing$Precision$estimate_values

valid_name <- config$wildcard_constraints$sample_id
if (!str_detect(sample_id, valid_name)) {stop('Sample id does not match supplied or default wildcard constraints!')}

use_chr <- get_chromosome_set()

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

############
# Run code #
############

pennCNVfiles <- c(paste0(sample_id, '/', sample_id, '.penncnv-auto.tsv'),
				  paste0(sample_id, '/', sample_id, '.penncnv-chrx.tsv'))
if (sex == 'm') { pennCNVfiles <- c(pennCNVfiles, paste0(sample_id, '/', sample_id, '.penncnv-chry.tsv')) }

results <- list()

if (args$penncnv) {
	res <- file.path(datapath, pennCNVfiles) %>%
		lapply(read_PennCNV) %>%
		bind_rows()
	if (nrow(res) > 0) {
	  results[['PennCNV']] <- res
	} else {
	  tool.order <- tool.order[tool.order != 'PennCNV']
	}
}
if (args$cbs) {
	res <- file.path(datapath, sample_id, paste0(sample_id, '.CBS.tsv')) %>%
		lapply(read_tsv, show_col_types = F) %>%
		bind_rows()
	if (nrow(res) > 0) {
	  results[['CBS']] <- res
	} else {
	  tool.order <- tool.order[tool.order != 'CBS']
	}
}


if (merge.before.filter) {
	results <- lapply(results, \(x) merge_calls(x, merge.distance)) %>%
		lapply(\(x) prefilter_calls(x, min.snp, min.length, min.snp.density))
} else {
	results <- lapply(results, prefilter_calls(x, min.snp, min.length, min.snp.density)) %>%
		lapply(\(x) merge_calls(x, merge.distance))
}


# check which calls overlap between tools, calculate coverage & convert cols to list, so they can store values from each tool post-overlap
# > for plotting original calls are still needed; these are retained (overlap.state == 'pre-overlap')
tools <- names(results)
if (!all(tool.order %in% tools)) {
	quit('"tool.order" contains tools that are not available')
} else if (!all(tools %in% tool.order)) {
	warning(paste0('Not all available CNV calling tools are specified in "tool.order".\n',
								 'These tools will *not* be used: ', paste0(tools[!tools %in% tool.order], collapse = ', ')))
}

processed.snps.file <- file.path(datapath,
                                sample_id,
								str_glue("{sample_id}.filtered-data-{n_snp.filterset}.tsv")
)
combined_tools_sample <- gr <- results[tools] %>%
	bind_ranges() %>%
	combine_CNV_callers(tool.overlap.greatest.call.min.perc, tool.overlap.median.cov.perc) %>%
	get_accurate_snp_probe_count(processed.snps.file)
 
if (!is.na(ref_id)) {
	
	fname <- file.path(datapath, ref_id, paste0(ref_id, '.combined-cnv-calls.tsv'))
	combined_tools_ref <- load_preprocessed_cnvs(fname) %>%
		# These cols will interefere
		dplyr::select(-length, -reference_overlap, -reference_coverage, -reference_caller,
									-n_genes, -overlapping_genes) %>%
		filter(caller_merging_state != 'pre-overlap') %>%
		as_granges()
	
	cnvs <- annotate_reference_overlap(combined_tools_sample, combined_tools_ref, min.reciprocal.coverage.with.ref)
	
} else {
	cnvs <- combined_tools_sample
}

cnvs <- cnvs %>%
	annotate_impact_lists(high_impact_gr, 'high_impact') %>%
	annotate_impact_lists(highlight_gr, 'highlight') %>%
	annotate_roi(sample_id, sampletable, gr_genes, gr_info) %>%
	annotate_gaps(config$static_data$array_gaps, min.perc.gap_area, gap_area.uniq_probes.rel) %>%
	annotate_high_density(config$static_data$array_density, density.quantile.cutoff) %>%
	finalise_gr_to_tb(gr_genes) %>%
	annotate_cnv.check.score(high_impact_gr, highlight_gr, check_scores) %>%
	annotate_precision.estimates(size_categories, precision_estimates)

		
outname.tsv <- file.path(datapath, sample_id, paste0(sample_id, '.combined-cnv-calls.tsv'))
cnvs.tb <- cnvs %>%
	rowwise() %>%
	mutate(across(one_of(list_cols), ~paste(., collapse=';')))
write_tsv(cnvs.tb, outname.tsv)
