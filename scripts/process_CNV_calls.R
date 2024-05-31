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
score_thresholds <- config$settings$CNV_processing$scoring_thresholds
impact_scores <- config$settings$CNV_processing$Impact_scoring
precision_extimates<- config$settings$CNV_processing$Precision_estimates

valid_name <- config$wildcard_constraints$sample_id
if (!str_detect(sample_id, valid_name)) {stop('Sample id does not match supplied or default wildcard constraints!')}

use_chr <- get_chromosome_set()

gr_genes <- load_gtf_data(config)
gr_info  <- load_genomeInfo(config)

high_impact_tb <- config$settings$CNV_processing$gene_overlap$high_impact_list %>%
	str_replace('__inbuilt__', config$snakedir) %>%
	read_tsv() %>%
	parse_hotspot_list()
highlight_tb <- config$settings$CNV_processing$gene_overlap$highlight_list %>%
	str_replace('__inbuilt__', config$snakedir) %>%
	read_tsv() %>%
	parse_hotspot_list()

########################
# Function definitions #
########################

## Merge
merge_calls <- function(df.or.GR) {
	message('Merging nearby raw calls')
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
	message('Pre-filtering calls')
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

# Caller overlap
combine_tools <- function(tools, min.greatest.region.overlap = 50, min.median.tool.coverage = 60) {
	message('Combining calls from multiple callers')
	
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
			              caller_merging_state = ifelse(plyranges::n() > 1, 'combined', 'no-overlap')) %>%
			as_tibble() %>%
			mutate(group.ID = paste('combined', CNV_type, seqnames, start, end, sep='_')) %>%
			dplyr::select(-strand, -seqnames) %>%
			rename_with(~paste0('group.', .), 1:3) %>%
			unnest(ID) %>%
			merge(as_tibble(gr), by = c('sample_id', 'CNV_type', 'ID')) %>%
			group_by(sample_id, CNV_type, pick(starts_with('group'))) %>%
			# The overlap here is *not* reciprocal overlap, as we need to deal with the possibility of one tool
			# 'merging' calls from another (or having >2 tools)
			mutate(overlap_merged_call = ifelse(caller_merging_state == 'combined', width / group.width * 100, NA),
			       max.ov = max(overlap_merged_call)) %>%
			group_by(sample_id, CNV_type, pick(starts_with('group')), CNV_caller) %>%
			# sum doesn't work if we have the filters in there and only group by CNV_caller
			mutate(tool.cov.sum = sum(overlap_merged_call)) %>%
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
			filter(caller_merging_state == 'combined' &
					   max.ov >= min.greatest.region.overlap &
					   tool.cov.ov.median >= min.median.tool.coverage
			) %>%
			select(-start, -end, -width,-ID) %>%
			rename_with(~str_remove(., '^group\\.')) %>%
			summarise(across(any_of(list_cols), ~list(.)),
					  overlap_merged_call = NA,
					  caller_merging_coverage = map2(CNV_caller, list(tool.cov.sum),
													\(CNV_caller, csum) tibble(t = CNV_caller, csum = csum) %>% unique() %>%
																	mutate(desc = paste0(t,'-', round(csum, 2))) %>%
																	pull(desc) %>% paste(collapse = ',')
													) %>% unlist(),
					  caller_merging_state = 'combined'
			) %>%
			as_granges() %>%
			ensure_list_cols()

		gr.unchanged <- ov %>%
			filter(caller_merging_state == 'no-overlap' |
					   max.ov < min.greatest.region.overlap |
					   tool.cov.ov.median < min.median.tool.coverage
			) %>%
			ungroup() %>%
			select(-starts_with('group.'), -max.ov, -tool.cov.sum, -tool.cov.ov.median) %>%
			mutate(caller_merging_coverage = NA_character_,
				   overlap_merged_call = NA,
			       caller_merging_state = 'no-overlap') %>%
			as_granges() %>%
			ensure_list_cols()

		gr.pre.overlap <-ov %>%
			filter(caller_merging_state == 'combined' &
					   max.ov >= min.greatest.region.overlap &
					   tool.cov.ov.median >= min.median.tool.coverage
			) %>%
			ungroup() %>%
			select(-starts_with('group.'), -max.ov, -tool.cov.sum, -tool.cov.ov.median) %>%
			mutate(caller_merging_coverage = NA_character_,
				   caller_merging_state = 'pre-overlap') %>%
			as_granges() %>%
			ensure_list_cols()

		return(bind_ranges(gr.changed, gr.unchanged, gr.pre.overlap))

	} else {
		gr <- gr %>%
			mutate(caller_merging_state = 'no-overlap',
				   #widths = NA,
				   overlap_merged_call = NA,
				   caller_merging_coverage = NA) %>%
			ensure_list_cols()

		return(gr)
	}

}

# n_nsp & uniq_snp annotation
get_accurate_snp_probe_count <- function(gr) {
	message('(re)calculating number of SNP probes per call')

	snp_probe_positions <- snp_probes_gr %>% reduce_ranges()

	n_snp <- count_overlaps(gr, snp_probes_gr)
	n_snp_uniq <- count_overlaps(gr, snp_probe_positions)

	gr$n_snp_probes <- n_snp
	gr$n_uniq_probe_positions <- n_snp_uniq

	gr

}


# Annotate (reciprocal'ish) overlap with calls from reference sample
annotate_ref_overlap <- function(gr_in, gr_ref, min.reciprocal.coverage.with.ref = 0.8) {
	message('Annotation calls with reference overlap')
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
			group_by(granges.x.seqnames, granges.x.start, granges.x.end, caller_merging_state.x, CNV_caller.x, CNV_type.x) %>%
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
			mutate(reference_overlap = FALSE,
				   reference_coverage = list(NA),
				   reference_caller = list(NA))
	}
	
	gr_out
	
}

# Annotate calls with hotspot / gene lists
annotate_impact_lists <- function(gr, list_name = 'high_impact') {
	message('Annotation calls with gene lists')

	if (list_name == 'high_impact') {
		hotspot_gr <- high_impact_tb
	} else if (list_name == 'highlight') {
		hotspot_gr <- highlight_tb
	} else {
		quit('Can only annotate "high_impact" or "highlight" lists')
	}

	# Make an extra col listing all overlapping directly defined genes
	gr@elementMetadata[[paste0(list_name, '_hits')]] <- NA_character_
	ov <- group_by_overlaps(gr, hotspot_gr)
	if (length(ov) > 0) {
		ov_hits <- ov %>%
			mutate(type_check = str_detect(CNV_type, str_replace(call_type, '^any$', '.*'))) %>%
			filter(type_check) %>%
			reduce_ranges(hits = paste(unique(hotspot),collapse = ','))
		gr[ov_hits$query,]@elementMetadata[[paste0(list_name, '_hits')]] <- ov_hits$hits
	}

	return(gr)
}


annotate_gaps <- function(gr, gapfile) {
	message('Annotation calls with gaps')

	gap_areas <- read_bed(gapfile) %>%
		select(-name, -score) %>%
		mutate(gap_size = width) %>%
		filter(seqnames %in% get_chromosome_set())

	# This assumes that calls _fully_ overlap gaps
	# which should be given since gaps are defined larte intre-probe distances
	gr$percent_gap_coverage <- join_overlap_left(gr, gap_areas) %>%
		group_by(ID) %>%
		reduce_ranges(gap_size_sum = sum(gap_size)) %>%
		mutate(perc_gap = ifelse(is.na(gap_size_sum), 0, gap_size_sum  / width)) %>%
		as_tibble() %>%
		pull(perc_gap)

	gap_slope <- gap_area.uniq_probes.rel[[1]]
	gap_intercept <- gap_area.uniq_probes.rel[[2]]

	gr <- gr %>%
		mutate(probe_coverage_gap =  ifelse(is.na(percent_gap_coverage),
										 FALSE,
										 percent_gap_coverage > min.perc.gap_area &
												(gap_slope * percent_gap_coverage + gap_intercept) <= log2(n_uniq_probe_positions))
		)

	return(gr)
}

annotate_high_density <- function(gr, density_file) {
	message('Annotation calls for high probe density')

	array_density <- read_bed(density_file) %>%
		select(-name) %>%
		mutate(density = score) %>%
		filter(seqnames %in% get_chromosome_set()) %>%
		filter(density > 0)

	density_quantile_cutoff <- quantile(array_density$density, density.quantile.cutoff)

	call_densities <- join_overlap_left(gr, array_density) %>%
		group_by(ID) %>%
		# partially overlap density windows might overproportionally affect calls this way
		# This only matters for very small calls thogub (and density in neighbouring windows will be similarish)
		reduce_ranges(density = mean(density)) %>%
		as_tibble() %>%
		pull(density)

	gr$high_probe_density <- call_densities > density_quantile_cutoff

	gr
}

annotate_roi <- function(gr, sample_id, sampletable) {

	if (! 'Regions_of_Interest' %in% colnames(sampletable)) {
		return(gr)
	}

	message('Annotating calls with ROI')

	roi <- sampletable[sampletable$Sample_ID == sample_id, ]$Regions_of_Interest

	if (is.na(roi) | is.null(roi) | roi == '') {
		return(gr)
	}

	roi_gr <- roi %>% str_split(';') %>% unlist() %>%
		as_tibble() %>% rename(roi = value) %>%
		mutate(hotspot = str_remove(roi, '^.*\\|'),
			   mapping = case_when(
				   str_detect(roi, '(chr)?[0-9XY]{1,2}:[0-9.,]+-[0-9.,]+') ~ 'position',
				   str_detect(roi, '[0-9XY]{1,2}(p|q)[0-9.]+') ~ 'gband',
				   TRUE ~ 'gene_name'
			   ),
			   roi_name = str_extract('^[^|]+\\|') %>% str_remove('\\|'),
			   roi_name = case_when(
				   mapping == 'gene_name' ~ hotspot,
				   is.na(roi_name) ~ paste0('ROI_', seq_along(roi_name)),
				   TRUE  ~ roi_name
			   ),
		) %>%
		parse_highligh_list()


	gr$ROI_hits <- NA_character_
	ov <- group_by_overlaps(gr, roi_gr)
	if (length(ov) > 0) {
		ov_hits <- ov %>%
			reduce_ranges(hits = paste(unique(roi_name),collapse = ','))
		gr[ov_hits$query,]$ROI_hits <- ov_hits$hits
	}

	gr
}


# Sanitize output & add gene overlap annotation
finalise_gr_to_tb <- function(gr) {

	gr$n_genes <- count_overlaps(gr, gr_genes)
	gr$overlapping_genes <- NA_character_
	ov_genes <- group_by_overlaps(gr, gr_genes) %>% reduce_ranges(genes = paste(gene_name, collapse = ','))
	gr[ov_genes$query,]$overlapping_genes <- ov_genes$genes

	tb <- as_tibble(gr) %>%
		rowwise() %>%
		mutate(across(any_of(list_cols), ~ list(.))) %>%
		bind_rows(expected_final_tb) %>%
		rowwise() %>%
		mutate(length = ifelse('lenght' %in% names(.), length, width)) %>%
		dplyr::select(one_of(colnames(expected_final_tb)))

	return(tb)
}

add_call_scoring <- function(tb) {

	tb %>%
	  rowwise() %>%
	  mutate(
		  size_category = case_when(
			length >= score_thresholds$extreme.loh & CNV_type %!in% c('gain', 'loss') ~ 'extreme',
			length >= score_thresholds$extreme.cnv & CNV_type %in% c('gain', 'loss') ~ 'extreme',
			length >= score_thresholds$very.large.loh & CNV_type %!in% c('gain', 'loss') ~ 'very_large',
			length >= score_thresholds$very.large.cnv & CNV_type %in% c('gain', 'loss') ~ 'very_large',
			length >= score_thresholds$large.loh & CNV_type %!in% c('gain', 'loss') ~ 'large',
			length >= score_thresholds$large.cnv & CNV_type %in% c('gain', 'loss') ~ 'large',
			length >= score_thresholds$medium.loh & CNV_type %!in% c('gain', 'loss') ~ 'medium',
			length >= score_thresholds$medium.cnv & CNV_type %in% c('gain', 'loss') ~ 'medium',
			TRUE ~ 'small'
		  ),
		  Impact_Score =
			ifelse(CNV_type %in% c('gain', 'loss'),
				   1/3 * log(length) * log(length) - 15,
				   0.275 * log(length) * log(length) - 15
			) +
			# Base score for hitting HI / HL / ROI
			(impact_scores$highimpact_base * !is.na(high_impact_hits) ) +
		    (impact_scores$highlight_base * !is.na(highlight_hits) ) +
		    (impact_scores$roi_hit_base * !is.na(ROI_hits) ) +
		    # Per gene/region score:
			# - scores per ROI (gene & others handeled the same)
			ifelse(!is.na(ROI_hits), (1 + str_count(ROI_hits, ',')) * impact_scores$pere_gene_roi, 0) +
			# - scores per non-gene HI / HL
			(high_impact_tb %>% as_tibble() %>%
				filter(mapping != 'gene_name' & hotspot %in% unlist(str_split(high_impact_hits, ','))) %>%
				mutate(impact_score = ifelse(is.na(impact_score), impact_scores$per_gene_highimpact, impact_score)) %>%
				pull(impact_score) %>%	sum()) +
			(highlight_tb %>% as_tibble() %>%
				filter(mapping != 'gene_name' & hotspot %in% unlist(str_split(highlight_hits, ','))) %>%
				mutate(impact_score = ifelse(is.na(impact_score), impact_scores$per_gene_highlight, impact_score)) %>%
				pull(impact_score) %>%	sum()) +
			# - score all genes that aren't also an ROI, use scores from tsv where available
			# dplyr will generate a bunch of warnings for calls without any genes
			suppressWarnings(bind_rows(
				high_impact_tb %>% as_tibble() %>%
					mutate(impact_score = ifelse(is.na(impact_score), impact_scores$per_gene_highimpact, impact_score)),
				highlight_tb %>% as_tibble() %>%
					mutate(impact_score = ifelse(is.na(impact_score), impact_scores$per_gene_highlight, impact_score)),
				str_split(overlapping_genes, ',') %>% unlist() %>%
					as_tibble() %>% rename(hotspot = value) %>%
					mutate(mapping = 'gene_name', impact_score = impact_scores$per_gene_any)
			  	) %>%
				filter(hotspot %in% unlist(str_split(overlapping_genes, ',')) &
						   hotspot %!in% unlist(str_split(ROI_hits, ',')) &
						   mapping == 'gene_name') %>%
				group_by(hotspot) %>% summarise(impact_score = max(impact_score)) %>%
				pull(impact_score) %>% sum())
		  ,
		  Precision_Estimate =
			  ifelse(CNV_type %in% c('gain', 'loss'),
				precision_extimates[[
					ifelse(caller_merging_state == 'combined', 'multiple_Callers', unlist(CNV_caller))]][[
					size_category]] +
				  (precision_extimates$Call_has_Gap * (probe_coverage_gap)) +
				  (precision_extimates$HighSNPDensity * (high_probe_density)),
				NA_real_
			  )
	  ) %>%
	  ungroup() %>%
	  dplyr::select(-size_category)

}


############
# Run code #
############

procssed.snps.file <- file.path(datapath,
                                sample_id,
								str_glue("{sample_id}.filtered-data-{n_snp.filterset}.tsv")
)
snp_probes_gr <- read_raw(procssed.snps.file) %>%
	as_granges(seqnames = Chr, start = Position, width = 1)

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

combined_tools_sample <- combine_tools(tools, tool.overlap.greatest.call.min.perc, tool.overlap.median.cov.perc) %>%
	get_accurate_snp_probe_count()
 
if (!is.na(ref_id)) {
	
	fname <- file.path(datapath, ref_id, paste0(ref_id, '.combined-cnv-calls.tsv'))
	combined_tools_ref <- load_preprocessed_cnvs(fname) %>%
		# These cols will interefere
		dplyr::select(-length, -reference_overlap, -reference_coverage, -reference_caller,
									-n_genes, -overlapping_genes) %>%
		filter(caller_merging_state != 'pre-overlap') %>%
		as_granges()
	
	cnvs <- annotate_ref_overlap(combined_tools_sample, combined_tools_ref, min.reciprocal.coverage.with.ref)
	
} else {
	cnvs <- combined_tools_sample
}

cnvs <- cnvs %>%
	annotate_impact_lists('high_impact') %>%
	annotate_impact_lists('highlight') %>%
	annotate_roi(sample_id, sampletable) %>%
	annotate_gaps(config$static_data$array_gaps) %>%
	annotate_high_density(config$static_data$array_density) %>%
	finalise_gr_to_tb() %>%
	add_call_scoring()

		
outname.tsv <- file.path(datapath, sample_id, paste0(sample_id, '.combined-cnv-calls.tsv'))
cnvs.tb <- cnvs %>%
	rowwise() %>%
	mutate(across(one_of(list_cols), ~paste(., collapse=';')))
write_tsv(cnvs.tb, outname.tsv)
