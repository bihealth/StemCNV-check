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

########################
# Function definitions #
########################

## Merge

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
	if (is.data.frame(df.or.GR)){
		df.or.GR <- dplyr::filter(df.or.GR, numsnp >= min.snp & length >= min.length & snp.density >= min.snp.density)
	} else {
		df.or.GR <- plyranges::filter(df.or.GR,	numsnp >= min.snp & width >= min.length & (numsnp / width * 1e6) >= min.snp.density )
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
			mutate(group.ID = paste('combined', CNV.state, seqnames, start, end, sep='_')) %>%
			dplyr::select(-strand, -seqnames) %>%
			rename_with(~paste0('group.', .), 1:3) %>%
			unnest(ID) %>%
			merge(as_tibble(gr), by = c('sample_id', 'CNV.state', 'ID')) %>%
			group_by(sample_id, CNV.state, pick(starts_with('group'))) %>%
			# The overlap here is *not* reciprocal overlap, as we need to deal with the possibility of one tool
			# 'merging' calls from another (or having >2 tools)
			mutate(coverage.overlap = ifelse(tool.overlap.state == 'combined', width / group.width * 100, NA),
			       max.ov = max(coverage.overlap)) %>%
			group_by(sample_id, CNV.state, pick(starts_with('group')), tool) %>%
			# sum doesn't work if we have the filters in there and only group by tool
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


# #################
#
# ##### Add VCF output #######
#
# #TODO add standard vcf output here instead of it's own script?
# # -> can run another script to get different vcf, maybe?
#
#
# # BCF2 is a binary, compressed equivalent of VCF that can be indexed
# # with tabix and can be efficiently decoded from disk or streams. For efficiency reasons BCF2 only supports a subset
# # of VCF, in that all info and genotype fields must have their full types specified
#
# # Could use cancer strategy to encode ref calls in vcf:
# ##PEDIGREE=<Derived=ID2,Original=ID1>
#
#
#
# # necessary info
# # ANY easy way to rebuild this?
#
#
# # Header content
#
# vcf.header <- c(
# 	'##fileformat=VCFv4.2',
# 	paste0('##fileDate=', Sys.Date()),
# 	#TODO: replace with proper name
# 	paste0('##source=', 'CNV-calling-pipeline'),
# 	paste0('##illumina-array-manifest=', str_remove(config$static_data$bpm_manifest_file, '\\.bpm')),
# 	'##ALT=<ID=DEL,Description="Deletion">',
# 	'##ALT=<ID=DUP,Description="Duplication">',
# 	'##ALT=<ID=LOH,Description="Loss of heterozygosity">'
# )
#
# #TODO: add contig info (i.e. from UCSC)
#
# chrom_levels <- get_chromosome_set(config)
# #TODO: get_chromosome_data()
#
# #need to check which chr format is used (both are allowed in principle
# ##contig=<ID=1,length=249250621,URL=...>
#
# info.lines <- c(
# 	'##INFO=<ID=END,Number=1,Type=Integer,Description="End coordinate of the variant">',
# 	'##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">',
# 	'##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">'
# )
#
# #TODO: think, does it make sense to include also the pre-filtered calls here?
# #TODO 2: this will be expanded on in the future; might become it's onw script
# filter.lines <- c(
# 	'##FILTER=<ID=PASS,Description="All filters passed">',
# 	'##FILTER=<ID=nonreportable,Description="CNV call does not match repotable QC treshholds">'
# )
#
#
# #? what about ref annotation ?
#
# #-> CNV/SV entries can do without proper original base info
#
# # Note that for deletions the position given is actually the base preceding the event.
#
# # POS is 1-indexed
#
# # REF can be N
# #  If any of the ALT alleles is a symbolic allele (an angle-bracketed ID String “<ID>”)
# #  then the padding base is required and POS denotes the coordinate of the base preceding the polymorphism
# #  the REF and ALT Strings must include the base before the event (which must be reflected in the POS field)
# #  --> can just use the 0-indexed numbers?
#
# # QUAL should be . with array data ( can't easily do a phred like score here )
#
# # INFO:
# # END is reserved & should be used; For precise variants, END is POS + length of REF allele - 1, and the for imprecise variants the corresponding best estimate.
#
# cnvs %>% 	mutate(
# 		seqnames = apply_chrom_format(), #str_remove(seqnames, remove_str),
# 		## Need to reduce listcols to a single value for vcf
# 		# Use highest tool_confidence  of single call
# 		call_confidence_score = get_a_combined_score_somehow(),
# 		# Max is not actually correct, but using sum here will be less accurate, this is not solvalble with the currently avalable information
# 		# Accurate number would need to be derived parallel to merging, requires loading the SNP subsets & doing a set merge on them
# 		#TODO: use the add_probe_count() function from validation analysis
# 		numsnp = numsnp %>% str_split(';') %>% unlist() %>% as.integer() %>% max(),
# 		tool = tool %>% str_split(';') %>% unlist() %>% unique() %>% sort %>% paste(collapse='&'),
# 		# No good way to decied if multiple copynumbers are predicted
# 		# Only difference here will by 3/4/... or 0/1
# 		# -> should maybe take the PennCNV one (instead of majority / 3 or 0)?
# 		# -> Maybe 3 or 0 are better assumptions if both 3/4 and 0/1 are the options.
# 		copynumber = copynumber %>% str_split(';') %>% unlist() %>%table() %>% which.max %>% names(),
# 		CNV.type = ifelse(CNV.state == 'gain', 'DUP', 'DEL'),
# 		CNV.type = ifelse(CNV.state == 'LOH', 'LOH', CNV.type),
# 		CHROM = factor(seqnames, levels = chrom_levels),
# 		POS = start, #TODO not sure if this is 0- or 1-indexed, check & adjist (should be 1-index for the base before the CNV)
# 		ID = str_glue("{CNV.type}_{seqnames}_{start}_{end}"),
# 		REF = 'N',
# 		ALT = str_glue('<{CNV.type}>'),
# 	    QUAL = '.',
# 		#TODO: thresholds are currently defined on report level, not globally; need something to filter here?
# 		# -> or change thresholds?
# 		FILTER = 'PASS',
# 		INFO = str_glue("END={end};SVTYPE={CNV.type};SVLEN={length}"),
# 		#needed?
# 		#FORMAT = 'CN:PN:TO:CO',
# )
#
#
#
# #TODO make it possible to configure this / add things in the future?
# vcf.format.content <- list(
# 	#Probably best to NOT include genotype -> possibly required for some use cases?
# 	#c('genotype', 'GT', "Segment genotype"),
# 	c('copynumber', 'CN', "Segment most-likely or estimated copy-number call"),
# 	c('numsnp', 'PN', "Number of SNP probes in the segment"),
# 	c('tool', 'TO', "CNV calling tools for segment"),
# 	#-> largely useless right now, but maybe we can find/define a better score?
# 	c('call_coonfidence_score ', 'CO', "Confidence score segment call (if available)")
# )
# vcf_cols <- c('CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT')
#
# write_to_vcf <- function(tb, outvcf) {
#
# 	format_str <- paste(sapply(vcf.format.content, function (x) x[2]), collapse = ':')
# 	use.cols <- sapply(vcf.format.content, function (x) x[1])
#
# 	format.lines <- sapply(vcf.format.content, function (x) {
# 		#Integer vs numeric ?!
# 		type <- ifelse(is.integer(tb[[x[1]]]), 'Integer', 'String')
# 		shorthand <- x[2]
# 		desc <- x[3]
# 		str_glue('##FORMAT=<ID={shorthand},NUMBER=1,Type={type},Description="{desc}">')
# 	})
#
# 	writeLines(c(vcf.header, info.lines, format.lines,
# 	             paste0('#', paste(c(vcf_cols, sample_id), collapse = '\t'))), con = outvcf)
#
# 	tb <- tb %>%
# 		filter(CNV.state %in% args$include_state)
#
# 	if (nrow(tb) > 0) {
# 		tb %>%
# 			rowwise() %>%
# 			mutate(
# 				FORMAT = format_str,
# 			) %>%
# 			do({
# 				out = as.data.frame(.)
# 				out$sample_formatted = out[,use.cols] %>% paste(collapse = ':')
# 				out
# 			}) %>%
# 			dplyr::select(one_of(c(vcf_cols, 'sample_formatted', 'sample_id'))) %>%
# 			pivot_wider(names_from = sample_id, values_from = sample_formatted) %>%
# 			arrange(CHROM, POS) %>%
# 			write_tsv(outvcf, append=T)
# 	}
#
# }
#
# #TODO: maybe should be an full cmd-line defined path?
# outvcf <- str_glue('{data_path}/{sample_id}/{sample_id}.combined-cnv-calls{name_addition}.vcf')
# processed.calls %>%
# 	filter(tool.overlap.state != 'pre-overlap') %>%
# 	#some purrr walk function instead?
# 	write_to_vcf(., outvcf = outvcf)
#
# #TODO
# # sys.call( gzip outfile )
# # sys.call( tabix outfile )