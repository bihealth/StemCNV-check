#! /usr/bin/Rscript
suppressMessages(library(argparse))

parser <- ArgumentParser(description="Collect data from CNV-calling tools and merge them into a single output, with added overlap stats")

parser$add_argument('data_path', type = 'character', help='Path to pipeline output data')
parser$add_argument('sample_id', type = 'character', help='Sample_ID to use')
parser$add_argument('configfile', type = 'character', help='Path to config file')
parser$add_argument('sampletablefile', type = 'character', help='Path to sampletable file')

parser$add_argument('-p', '--penncnv', action="store_true", help="Use PennCNV calls")
parser$add_argument('-c', '--cbs', action="store_true", help="Use CBS calls")
parser$add_argument('-g', '--gada', action="store_true", help="Use GADA calls")

# args <- parser$parse_args(c("-p", "-c", "test/data", "BIHi250-A-2", "test/test_config.yaml", "sample_table_example.txt"))
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

########################
# Function definitions #
########################

## Merge

#Wishlist: try this instead of distance based merge
# get_distance_to_probe <- function(position, n_probes) {
# 	1
# }

merge_calls <- function(df.or.GR) {
	if (is.data.frame(df.or.GR)) {
		df.or.GR <- as_granges(df.or.GR, seqnames = Chr)
	}
	df.or.GR %>%
		group_by(CNV.state, sample_id, tool) %>%
		stretch(merge.distance) %>%
		reduce_ranges(
			merged_tool_calls = plyranges::n(),
			numsnp = sum(numsnp),
			copynumber = paste(unique(copynumber), collapse = ','),
			tool_confidence = median(tool_confidence)
			) %>%
		stretch(-1*merge.distance) %>%
		mutate(ID = paste(tool, CNV.state, seqnames, start, end, sep='_'))
}

## Pre-filter
prefilter_calls <- function(df.or.GR) {
	# 'Our' chip has ~2.5 snps per 10kb -> average expected 250 per Mb
    # Ideally a threshold on desnity would apply to a pre-calculated window based SNP density
	# -> this could be done for each array & included in report but also used for i.e. PennCNV extra data (like GC model)
	if (is.data.frame(df.or.GR)){
		df.or.GR <- dplyr::filter(df.or.GR, numsnp >= min.snp & length >= min.length & snp.density >= min.snp.density) # &
										#(numsnp >= min.snp.or.length | length >= min.length.or.snp))
	} else {
		df.or.GR <- plyranges::filter(df.or.GR,	numsnp >= min.snp & width >= min.length & (numsnp / width * 1e6) >= min.snp.density )#&
										#(numsnp >= min.snp.or.length | width >= min.length.or.snp))
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
		group_by(sample_id, CNV.state) %>%
		reduce_ranges()

	# If overlaps exist (if not reduce_ranges can't make proper list cols):
	# Merge any number of overlaps for the same CNV.state together, summarise metadata into list cols
	# & calculate coverage of the final merged call by individual ones	
	if (length(ov_test) < length(gr))  {
		ov <- gr %>%
			group_by(sample_id, CNV.state) %>%
			reduce_ranges(ID = ID,
			              tool.overlap.state = ifelse(plyranges::n() > 1, 'combined', 'no-overlap')) %>%
			as_tibble() %>%
			# TODO: add numsnp for new grouped size -> need to read filtered probes for that!
			mutate(group.ID = paste('combined', CNV.state, seqnames, start, end, sep='_')) %>%
			dplyr::select(-strand, -seqnames) %>%
			rename_with(~paste0('group.', .), 1:3) %>%
			unnest(ID) %>%
			merge(as_tibble(gr), by = c('sample_id', 'CNV.state', 'ID')) %>%
			group_by(sample_id, CNV.state, pick(starts_with('group'))) %>%
			mutate(coverage.overlap = ifelse(tool.overlap.state == 'combined', width / group.width * 100, NA),
			       max.ov = max(coverage.overlap)) %>%
			group_by(sample_id, CNV.state, pick(starts_with('group')), tool) %>%
			# sum doesn't work if we have the filters in there and only group by too
			mutate(tool.cov.sum = sum(coverage.overlap)) %>%
			group_by(sample_id, CNV.state, seqnames, pick(starts_with('group'))) %>%
			mutate(tool.cov.ov.median = map2(list(tool), list(tool.cov.sum),
												\(tool, csum) tibble(tool = tool, csum = csum) %>% unique() %>%
																pull(csum) %>% median()
												),
			)

		#Now filter based on overlap of tool calls with merged region to check if the overlap can be accepted
		# - require at least 50% of the merged region to be covered by a single call to prevent chained overlaps
		# - require that the combined regions coverage from tools has a median of at least 60% (avg for only 2 tools)
		# TODO (for more than 2 tools): add some criteria to remove individual tools/calls from an combined region
		#  to see if the combined group of two other tools can be 'rescued'
		gr.changed <- ov %>%
			filter(tool.overlap.state == 'combined' &
					   max.ov >= min.greatest.region.overlap &
					   tool.cov.ov.median >= min.median.tool.coverage
			) %>%
			select(-start, -end, -width,-ID) %>%
			rename_with(~str_remove(., '^group\\.')) %>%
			summarise(across(any_of(list_cols), ~list(.)),
					  #tool is already a list now!
					  tool.coverage.overlap = map2(tool, list(tool.cov.sum),
													\(tool, csum) tibble(t = tool, csum = csum) %>% unique() %>%
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
				   coverage.by.ref = round(100 * width.overlap / granges.x.width, 2),
				   cov.ref.by.sample = round(100 * width.overlap / granges.y.width, 2),
				   # need to flatten/fully unpack tool.y or we get list of lists
				   ref.tool = list(unlist(tool.y)),
				   ref.state = CNV.state.y
			) %>%
			dplyr::select(-contains('.y'), -contains('width')) %>%
			filter(coverage.by.ref >= min.reciprocal.coverage.with.ref &
					   cov.ref.by.sample >= min.reciprocal.coverage.with.ref &
					   CNV.state.x == ref.state)
	)
	
	if (nrow(gr) > 0) {
		gr <- gr %>%
			group_by(granges.x.seqnames, granges.x.start, granges.x.end, tool.overlap.state.x, tool.x, CNV.state.x) %>%
			summarise(across(ends_with('.x'), ~ unique(.)),
								call.in.reference = TRUE,
								coverage.by.ref = list(coverage.by.ref),
								ref.tool = list(unlist(ref.tool)),
			) %>%
			dplyr::rename_with(~ str_remove(., '.x$') %>% str_remove('^granges.x.')) %>%
			as_granges()

		# `gr` only has (filtered) overlaps, need to rebuild the full callset if matching ref calls were found
		# Get Original regions not in the merged call set by ID
		non.ovs <- gr_in %>% filter(ID %!in% gr$ID)
		# mutate fails on empty GRanges object
		if (length(non.ovs) > 0) { non.ovs <- non.ovs %>% mutate(call.in.reference = FALSE) }

		gr_out <- bind_ranges(
			# Calls with matching reference
			gr,
			# Remaining calls not (sufficiently) overlapping reference
			non.ovs
		) 
		
	} else {
		gr_out <- gr_in %>%
			mutate(call.in.reference = FALSE,
						 coverage.by.ref = list(NA),
						 ref.tool = list(NA))
	}
	
	gr_out
	
}

# Sanitize output & add gene overlap annotation
finalise_tb <- function(gr) {
	
	gr$n_genes <- count_overlaps(gr, gr_genes)
	gr$overlapping.genes <- NA_character_
	ov_genes <- group_by_overlaps(gr, gr_genes) %>% reduce_ranges(genes = paste(gene_name, collapse = ','))
	gr[ov_genes$query,]$overlapping.genes <- ov_genes$genes
	
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
if (args$gada) {
	res <- file.path(datapath, sample_id, paste0(sample_id, '.GADA.tsv')) %>%
		lapply(read_GADA) %>%
		bind_rows()
	if (nrow(res) > 0) {
	  results[['GADA']] <- res
	} else {
	  tool.order <- tool.order[tool.order != 'GADA']
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
		dplyr::select(-length, -call.in.reference, -coverage.by.ref, -ref.tool,
									-n_genes, -overlapping.genes) %>%
		filter(tool.overlap.state != 'pre-overlap') %>%
		as_granges()
	
	cnvs <- annotate_ref_overlap(combined_tools_sample, combined_tools_ref, min.reciprocal.coverage.with.ref) %>%
		finalise_tb()
	
} else {
	cnvs <- combined_tools_sample %>%
		finalise_tb()	
}
		
outname.tsv <- file.path(datapath, sample_id, paste0(sample_id, '.combined-cnv-calls.tsv'))
cnvs.tb <- cnvs %>%
	rowwise() %>%
	mutate(across(one_of(list_cols), ~paste(., collapse=';')))
write_tsv(cnvs.tb, outname.tsv)

#TODO add standard vcf output here instead of it's own script?
# -> can run another script to get different vcf, maybe?
