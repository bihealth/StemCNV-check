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
		return ('scripts')
	}
}
source(file.path(get_script_dir(), 'R_io_functions.R'))

##################
# Variable setup #
##################

datapath <- args$data_path
sample_id <- args$sample_id
config <- read_yaml(args$configfile)
sampletable <- read_tsv(args$sampletablefile, col_types = 'cccccc', comment = '#')

sex <- get_sample_info(sample_id, 'sex', sampletable)
ref_id <- get_sample_info(sample_id, 'ref_id', sampletable)
sex.ref <- get_sample_info(sample_id, 'sex.ref', sampletable)

# Call merging & pre-filtering
min.snp				<- config$settings$postprocessing$min.snp
min.length			<- config$settings$postprocessing$min.length
min.snp.density 	<- config$settings$postprocessing$min.snp.density
merge.distance  	<- config$settings$postprocessing$merge.distance
merge.before.filter <- config$settings$postprocessing$merge.before.filter
# Tool overlapping
tool.order <- config$settings$CNV.calling.tools
tool.overlap.greatest.call.min.perc <- config$settings$postprocessing$tool.overlap.greatest.call.min.perc
tool.overlap.median.cov.perc <- config$settings$postprocessing$tool.overlap.median.cov.perc
# Compare to reference
min.reciprocal.coverage.with.ref <- config$settings$postprocessing$min.reciprocal.coverage.with.ref

valid_name <- config$wildcard_constraints$sample_id
if (!str_detect(sample_id, valid_name)) {stop('Sample id does not match supplied or default wildcard constraints!')}

use_chr <- get_chromosome_set()

gr_genes <- load_gtf_data(config)
gr_info  <- load_genomeInfo(config)

########################
# Function definitions #
########################

## Merge

merge_calls <- function(df.or.GR) {
	if (is.data.frame(df.or.GR)) {
		df.or.GR <- as_granges(df.or.GR, seqnames = Chr)
	}
	df.or.GR %>%
		group_by(CNV_type, sample_id, CNV_caller) %>%
		stretch(merge.distance) %>%
		reduce_ranges(
			n_premerged_calls = plyranges::n(),
			n_snp_probes = sum(n_snp_probes),
			copynumber = paste(unique(copynumber), collapse = ','),
			caller_confidence = median(caller_confidence)
			) %>%
		stretch(-1*merge.distance) %>%
		mutate(ID = paste(CNV_caller, CNV_type, seqnames, start, end, sep='_'))
}

## Pre-filter
prefilter_calls <- function(df.or.GR) {
	if (is.data.frame(df.or.GR)){
		df.or.GR <- dplyr::filter(df.or.GR, n_snp_probes >= min.snp & length >= min.length & snp.density >= min.snp.density)
	} else {
		df.or.GR <- plyranges::filter(df.or.GR,	n_snp_probes >= min.snp & width >= min.length & (n_snp_probes / width * 1e6) >= min.snp.density )
	}
	df.or.GR
}

ensure_list_cols <- function(tb.or.gr){
	as_tibble(tb.or.gr) %>%
		rowwise() %>%
		mutate(across(any_of(list_cols), ~ list(.))) %>%
		as_granges()
}

combine_tools <- function(tools, min.greatest.region.overlap = 50, min.median.tool.coverage = 60) {
	
	gr <- results[tools] %>%
		bind_ranges()
	
	ov_test <- gr %>%
		group_by(sample_id, CNV_type) %>%
		reduce_ranges()

	# If overlaps exist (if not reduce_ranges can't make proper list cols):
	# Merge any number of overlaps for the same CNV.state together, summarise metadata into list cols
	# & calculate coverage of the final merged call by individual ones	
	if (length(ov_test) < length(gr))  {
		ov <- gr %>%
			group_by(sample_id, CNV_type) %>%
			reduce_ranges(ID = ID,
			              tool.overlap.state = ifelse(plyranges::n() > 1, 'combined', 'no-overlap')) %>%
			as_tibble() %>%
			mutate(group.ID = paste('combined', CNV_type, seqnames, start, end, sep='_')) %>%
			dplyr::select(-strand, -seqnames) %>%
			rename_with(~paste0('group.', .), 1:3) %>%
			unnest(ID) %>%
			merge(as_tibble(gr), by = c('sample_id', 'CNV_type', 'ID')) %>%
			group_by(sample_id, CNV_type, pick(starts_with('group'))) %>%
			# The overlap here is *not* reciprocal overlap, as we need to deal with the possibility of one tool
			# 'merging' calls from another (or having >2 tools)
			mutate(coverage.overlap = ifelse(tool.overlap.state == 'combined', width / group.width * 100, NA),
			       max.ov = max(coverage.overlap)) %>%
			group_by(sample_id, CNV_type, pick(starts_with('group')), CNV_caller) %>%
			# sum doesn't work if we have the filters in there and only group by CNV_caller
			mutate(tool.cov.sum = sum(coverage.overlap)) %>%
			group_by(sample_id, CNV_type, seqnames, pick(starts_with('group'))) %>%
			mutate(tool.cov.ov.median = map2(list(CNV_caller), list(tool.cov.sum),
												\(CNV_caller, csum) tibble(CNV_caller = CNV_caller, csum = csum) %>% unique() %>%
																pull(csum) %>% median()
												),
			)

		#Now filter based on overlap of tool(CNV_caller) calls with merged region to check if the overlap can be accepted
		# - require at least 50% of the merged region to be covered by a single call to prevent chained overlaps
		# - require that the combined regions coverage from tools has a median of at least 60% (avg for only 2 tools)
		gr.changed <- ov %>%
			filter(tool.overlap.state == 'combined' &
					   max.ov >= min.greatest.region.overlap &
					   tool.cov.ov.median >= min.median.tool.coverage
			) %>%
			select(-start, -end, -width,-ID) %>%
			rename_with(~str_remove(., '^group\\.')) %>%
			summarise(across(any_of(list_cols), ~list(.)),
					  #tool is already a list now!
					  tool.coverage.overlap = map2(CNV_caller, list(tool.cov.sum),
													\(CNV_caller, csum) tibble(t = CNV_caller, csum = csum) %>% unique() %>%
																	mutate(desc = paste0(t,'-', round(csum, 2))) %>%
																	pull(desc) %>% paste(collapse = ',')
													) %>% unlist(),
					  tool.overlap.state = 'combined'
			) %>%
			as_granges() %>%
			ensure_list_cols()

		gr.unchanged <- ov %>%
			filter(tool.overlap.state == 'no-overlap' |
					   max.ov < min.greatest.region.overlap |
					   tool.cov.ov.median < min.median.tool.coverage
			) %>%
			ungroup() %>%
			select(-starts_with('group.'), -max.ov, -tool.cov.sum, -tool.cov.ov.median) %>%
			mutate(tool.coverage.overlap = NA_character_,
			       tool.overlap.state = 'no-overlap') %>%
			as_granges() %>%
			ensure_list_cols()

		gr.pre.overlap <-ov %>%
			filter(tool.overlap.state == 'combined' &
					   max.ov >= min.greatest.region.overlap &
					   tool.cov.ov.median >= min.median.tool.coverage
			) %>%
			ungroup() %>%
			select(-starts_with('group.'), -max.ov, -tool.cov.sum, -tool.cov.ov.median) %>%
			mutate(tool.coverage.overlap = NA_character_,
				   overlap.coverage = list(NA),
				   tool.overlap.state = 'pre-overlap') %>%
			as_granges() %>%
			ensure_list_cols()

		return(bind_ranges(gr.changed, gr.unchanged, gr.pre.overlap))

	} else {
		gr <- gr %>%
			mutate(tool.overlap.state = 'no-overlap',
				   #widths = NA,
				   coverage.overlap = NA,
				   tool.coverage.overlap = NA) %>%
			ensure_list_cols()

		return(gr)
	}

}


annotate_ref_overlap <- function(gr_in, gr_ref, min.reciprocal.coverage.with.ref = 0.8) {

	#If pair_overlaps is an empty set/df the min/max functions will give a warning, while later Granges functions would crash
	# dplyr::summarise will also give a (deprecation) warning if it gets a fully empty tibble
	suppressWarnings(
	gr <- pair_overlaps(gr_in, gr_ref) %>%
			as_tibble() %>% rowwise() %>%
			mutate(width.combined = max(granges.x.end, granges.y.end) - min(granges.x.start, granges.y.start) + 1,
				   width.overlap = min(granges.x.end, granges.y.end) - max(granges.x.start, granges.y.start) + 1,
				   reference_coverage = round(100 * width.overlap / granges.x.width, 2),
				   cov.ref.by.sample = round(100 * width.overlap / granges.y.width, 2),
				   # need to flatten/fully unpack tool.y or we get list of lists
				   reference_caller = list(unlist(CNV_caller.y)),
				   ref.state = CNV_type.y
			) %>%
			dplyr::select(-contains('.y'), -contains('width')) %>%
			filter(reference_coverage >= min.reciprocal.coverage.with.ref &
					   cov.ref.by.sample >= min.reciprocal.coverage.with.ref &
					   CNV_type.x == ref.state)
	)
	
	if (nrow(gr) > 0) {
		gr <- gr %>%
			group_by(granges.x.seqnames, granges.x.start, granges.x.end, tool.overlap.state.x, CNV_caller.x, CNV_type.x) %>%
			summarise(across(ends_with('.x'), ~ unique(.)),
					  reference_overlap = TRUE,
					  reference_coverage = list(reference_coverage),
					  reference_caller = list(unlist(reference_caller)),
			) %>%
			dplyr::rename_with(~ str_remove(., '.x$') %>% str_remove('^granges.x.')) %>%
			as_granges()

		# `gr` only has (filtered) overlaps, need to rebuild the full callset if matching ref calls were found
		# Get Original regions not in the merged call set by ID
		non.ovs <- gr_in %>% filter(ID %!in% gr$ID)
		# mutate fails on empty GRanges object
		if (length(non.ovs) > 0) { non.ovs <- non.ovs %>% mutate(reference_overlap = FALSE) }

		gr_out <- bind_ranges(
			# Calls with matching reference
			gr,
			# Remaining calls not (sufficiently) overlapping reference
			non.ovs
		) 
		
	} else {
		gr_out <- gr_in %>%
			mutate(call.in.reference = FALSE,
				   reference_coverage = list(NA),
				   reference_caller = list(NA))
	}
	
	gr_out
	
}

parse_highligh_list <- function(yaml_obj) {

	if (!is.null(yaml_obj$gene_symbol)) {
		gr_g   <- gr_genes %>%
			mutate(gene_hit = TRUE) %>%
			filter(gene_name %in% yaml_obj$gene_symbol)
	} else {
		gr_g <- GRanges(gene_name = character(),
						gene_hit = logical())
	}

	if (!is.null(yaml_obj$position)) {
		regex <- paste0('^(',
						paste(str_replace(yaml_obj$position, fixed('.'), '\\.'), collapse='|'),
						')')
		gr_pos <- gr_info %>%
			mutate(gene_hit = FALSE) %>%
			filter(str_detect(section_name, regex))
	} else {
		gr_pos <- GRanges(section_name = character(),
		                  gene_hit = logical())
	}

	bind_ranges(gr_g, gr_pos)

}


annotate_impact_lists <- function(gr, config, lists = 'high_impact') {

	if (lists == 'high_impact') {
		config_sect <- config$settings$gene_overlap$highimpact_curated_lists
		col_name <- 'high_impact'
	} else if (lists == 'highlight') {
		config_sect <- config$settings$gene_overlap$highlight_lists
		col_name <- 'highlight'
	} else {
		quit('Can only annotate "high_impact" or "highlight" lists')
	}
	col_name_g <- paste0(col_name, '_genes')

	yaml_files <- unlist(config_sect) %>%
		str_replace('__inbuilt__', config$snakedir)

	# Make a GRList object with coordinates for each gene/position of the lists
	highlight_lists <- lapply(yaml_files, read_yaml)
	names(highlight_lists) <- sapply(highlight_lists, \(x) x$name)

	highlight_gr <- lapply(highlight_lists, parse_highligh_list) %>%
		GRangesList()
	# Make a list-col & add highlight-list names that have any hits
	gr@elementMetadata[[col_name]] <- imap(highlight_gr, function(x, name) ifelse(count_overlaps(gr, x)>0, name,NA_character_ )) %>%
		pmap(\(...) na.omit(c(...)) %>% as.character())

	# Make an extra col listing all overlapping directly defined genes
	gr@elementMetadata[[col_name_g]] <- NA_character_
	ov_hits <- group_by_overlaps(gr, unlist(highlight_gr)) %>%
		reduce_ranges(genes = paste(unique(gene_name[gene_hit]),collapse = ','),
		              pos = paste(unique(section_name[!gene_hit]),collapse = ',')) %>%
		mutate(hits = paste(genes, pos, sep=',') %>% str_remove(',$'))
	gr[ov_hits$query,]@elementMetadata[[col_name_g]] <- ov_hits$hits

	return(gr)
}


annotate_gaps <- function(gr, gapfile) {

	gaps <- read_bed(gapfile) %>%
		select(-name, -score) %>%
		mutate(gap_size = width) %>%
		filter(seqnames %in% get_chromosome_set())

	# TODO calculated gap as in paper analysis (% of calls size vs n_probes)
	# gap_designation =  ifelse(is.na(perc_gap), F, perc_gap > 1/3 & (gap_slope * perc_gap + gap_intercept) <= log2(region.n_positions)

	#TODO this needs proper n_probes
	join_overlap_left(gr, gaps) %>%
	   group_by(ID, n_probes) %>%
	   reduce_ranges(gap_size_sum = sum(gap_size)) %>%
	   mutate(perc_gap = ifelse(is.na(gap_size_sum), 0, gap_size_sum  / width)) %>%
	   pull()

	gr <- gr %>%


		#mutate(gap = ifelse(count_overlaps(gr, gaps) > 0, TRUE, FALSE))
	return(gr)
}


# Sanitize output & add gene overlap annotation
finalise_tb <- function(gr) {
	
	gr$n_genes <- count_overlaps(gr, gr_genes)
	gr$overlapping_genes <- NA_character_
	ov_genes <- group_by_overlaps(gr, gr_genes) %>% reduce_ranges(genes = paste(gene_name, collapse = ','))
	gr[ov_genes$query,]$overlapping_genes <- ov_genes$genes
	
	tb <- as_tibble(gr) %>%
		rowwise() %>%
		mutate(across(one_of(list_cols), ~ list(.))) %>%
		bind_rows(expected_final_tb) %>%
		rowwise() %>%
		mutate(length = ifelse('lenght' %in% names(.), length, width)) %>%
		dplyr::select(one_of(colnames(expected_final_tb))) 
	
	return(tb)
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
	results <- lapply(results, merge_calls) %>%
		lapply(prefilter_calls)
} else {
	results <- lapply(results, prefilter_calls) %>%
		lapply(merge_calls)
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

combined_tools_sample <- combine_tools(tools, tool.overlap.greatest.call.min.perc, tool.overlap.median.cov.perc)
 
if (!is.na(ref_id)) {
	
	fname <- file.path(datapath, ref_id, paste0(ref_id, '.combined-cnv-calls.tsv'))
	combined_tools_ref <- load_preprocessed_cnvs(fname) %>%
		# These cols will interefere
		dplyr::select(-length, -reference_overlap, -reference_coverage, -reference_caller,
									-n_genes, -overlapping_genes) %>%
		filter(tool.overlap.state != 'pre-overlap') %>%
		as_granges()
	
	cnvs <- annotate_ref_overlap(combined_tools_sample, combined_tools_ref, min.reciprocal.coverage.with.ref) %>%
		annotate_impact_lists(config, 'high_impact') %>%
		annotate_impact_lists(config, 'highlight') %>%
		finalise_tb()
	
} else {
	cnvs <- combined_tools_sample %>%
		annotate_impact_lists(config, 'high_impact') %>%
		annotate_impact_lists(config, 'highlight') %>%
		finalise_tb()
}
		
outname.tsv <- file.path(datapath, sample_id, paste0(sample_id, '.combined-cnv-calls.tsv'))
cnvs.tb <- cnvs %>%
	rowwise() %>%
	mutate(across(one_of(list_cols), ~paste(., collapse=';')))
write_tsv(cnvs.tb, outname.tsv)
