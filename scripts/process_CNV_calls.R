#! /usr/bin/Rscript
# Filter SNP probes by quality scores
suppressMessages(library(argparse))

parser <- ArgumentParser(description="Collect data from CNV-calling tools and merge them into a single output, with added overlap stats")

parser$add_argument('data_path', type = 'character', help='Path to pipeline output data')
parser$add_argument('sample_id', type = 'character', help='Sample_ID to use')
parser$add_argument('configfile', type = 'character', help='Path to config file')
parser$add_argument('sampletablefile', type = 'character', help='Path to sampletable file')

parser$add_argument('-p', '--penncnv', action="store_true", help="Use PennCNV calls")
parser$add_argument('-c', '--cbs', action="store_true", help="Use CBS calls")
parser$add_argument('-g', '--gada', action="store_true", help="Use GADA calls")

args <- parser$parse_args()
# args <- parser$parse_args(c("-p", "-c", "data_new", "BIHi001-B_WB02", "/tmp/tmpewl9ym4n.yaml", "sample_table.txt"))
# args <- parser$parse_args(c("-p","-c", "/data/cephfs-1/work/groups/cubi/projects/2023-04-24_Schaefer_Array_CNV/data", "9807821008_R02C01", "/data/gpfs-1/users/vonkunic_c/scratch/tmp/hpc-cpu-63/tmpxcb94bht.yaml", "sample_table.txt"))


suppressMessages(library(plyranges))
suppressMessages(library(GenomicRanges))
suppressMessages(library(tidyverse))
suppressMessages(library(yaml))


##################
# Variable setup #
##################

datapath <- args$data_path
sample_id <- args$sample_id
config <- read_yaml(args$configfile)
sampletable <- read_tsv(args$sampletablefile, col_types = 'cccccc', comment = '#')

sex <- sampletable[sampletable$Sample_ID == sample_id, ]$Sex %>%
	tolower() %>% substr(1,1)

ref_id <- sampletable[sampletable$Sample_ID == sample_id, ]$Reference_Sample

if (!is.na(ref_id)){
  sex.ref <- sampletable[sampletable$Sample_ID == ref_id, ]$Sex %>%
  	tolower() %>% substr(1,1)

  if(sex.ref != sex) {
  	stop('Sex of sample and reference does not match!')
  }
}

# Call merging & pre-filtering
min.snp				<- config$settings$postprocessing$min.snp
min.length			<- config$settings$postprocessing$min.length
min.snp.density 	<- config$settings$postprocessing$min.snp.density
min.snp.or.length	<- config$settings$postprocessing$min.snp.or.length
min.length.or.snp	<- config$settings$postprocessing$min.length.or.snp
merge.distance  	<- config$settings$postprocessing$merge.distance
merge.before.filter <- config$settings$postprocessing$merge.before.filter
# Tool overlapping
tool.order <- config$settings$CNV.calling.tools
tool.overlap.min.perc <- config$settings$postprocessing$tool.overlap.min.perc
# Compare to reference
min.cnv.coverage.by.ref <- config$settings$postprocessing$min.cnv.coverage.by.ref

valid_name <- config$wildcard_constraints$sample_id
# valid_name=ifelse(is.null(valid_name), '[0-9]{12}_R[0-9]{2}C[0-9]{2}', valid_name)
if (!str_detect(sample_id, valid_name)) {stop('Sample id does not match supplied or default wildcard constraints!')}

use_chr <- paste0('chr', c(1:22, 'X', 'Y'))

gtf_file <- config$static_data$genome_gtf_file
exclude_regexes <- config$settings$gene_overlap$exclude_gene_type_regex %>%
        paste(collapse = '|')
gene_type_whitelist <- config$settings$gene_overlap$include_only_these_gene_types

gr_genes  <- read_gff(gtf_file, col_names = c('source', 'type', 'gene_id', 'gene_type', 'gene_name')) %>%
    filter(type == 'gene' & seqnames %in% use_chr & !str_detect(gene_type, exclude_regexes))
if (typeof(gene_type_whitelist) == 'list' & length(gene_type_whitelist) > 0){
    gr_genes <- filter(gr_genes, gene_type %in% gene_type_whitelist)
}

########################
# Function definitions #
########################

## PennCNV
read_PennCNV <- function(filename) {
	read.table(filename, sep='', header = F, fill=T,
						 col.names = c('Position', 'numsnp', 'length', 'hmm.state', 'input', 'startsnp', 'endsnp', 'conf')) %>% 
		separate(Position, c('Chr', 'start_pos', 'end_pos'), convert=T) %>%
		dplyr::rename(start = start_pos, end = end_pos, sample_id = input) %>%
		mutate(across(c(4,5,8,9,10), ~ str_remove(., '.*=')),
					 across(c(4,10), ~as.numeric(.)),  
					 Chr = factor(Chr, levels = c(paste0('chr', 1:22), 'chrX', 'chrY')),
					 length = str_remove_all(length, ',') %>% as.integer(),
					 snp.density = numsnp / length * 1e6,
					 copynumber = str_extract(hmm.state, '(?<=cn=)[0-9]') %>% as.integer(),
					 hmm.state = str_remove(hmm.state, ',cn=[0-9]'),
					 CNV.state = ifelse(copynumber < 2, 'loss', NA),
					 CNV.state = ifelse(copynumber == 2, 'LOH', CNV.state),
					 CNV.state = ifelse(copynumber > 2, 'gain', CNV.state),
					 CNV.state = as.character(CNV.state),
					 tool = 'PennCNV',
			         # basename can't handly empty input ...
					 sample_id = str_remove(sample_id, '.*/') %>% str_remove('\\.filtered-data-.*\\.tsv$'),
		)
}

## CBS
read_CBS <- function(filename) {
	tb <- read_tsv(filename, show_col_types = F)
	tb
}


## GADA / MAD

read_GADA <- function(filename) {
	read_tsv(filename, show_col_types = F) %>%
		mutate(conf = NA,
			   snp.density = numsnp / length * 1e6,
			   copynumber = ifelse(!CNV.state %in% c('gain', 'loss'), 2, 3),
			   copynumber = ifelse(CNV.state == 'loss', 1, copynumber),)
}


## Merge

#TODO: try this instead of distance based merge
get_distance_to_probe <- function(position, n_probes) {
	1
	
}

merge_calls <- function(df.or.GR) {
	if (is.data.frame(df.or.GR)) {
		df.or.GR <- makeGRangesFromDataFrame(df.or.GR, keep.extra.columns = T, ignore.strand = T, seqnames.field = 'Chr', start.field = 'start', end.field = 'end')
	}
	df.or.GR %>%
		group_by(CNV.state, sample_id, tool) %>%
			stretch(merge.distance) %>%
			reduce_ranges(
				individual_calls = plyranges::n(),
				numsnp = sum(numsnp),
				copynumbers = paste(unique(copynumber), collapse = ','),
				conf = median(conf)
				) %>%
			stretch(-1*merge.distance)
}

## Pre-filter
prefilter_calls <- function(df.or.GR) {
	if (is.data.frame(df.or.GR)){
		df.or.GR <- dplyr::filter(df.or.GR, 
									numsnp >= min.snp & length >= min.length & snp.density >= min.snp.density &
										(numsnp >= min.snp.or.length | length >= min.length.or.snp))
	} else {
		df.or.GR <- plyranges::filter(df.or.GR,
									numsnp >= min.snp & width >= min.length & (numsnp / width * 1e6) >= min.snp.density &
										(numsnp >= min.snp.or.length | width >= min.length.or.snp))
	}
	df.or.GR
}

#TODO: move this + input function defs to another script?
# Harmonize output
expected_final_tb <- tibble(
	sample_id = character(),
	seqnames = factor(c(), levels = use_chr),
	start = integer(),
	end = integer(),
	length = integer(),
	CNV.state = character(),
	call.in.reference = logical(),
	coverage.by.ref = list(),
	tool = list(),
	individual_calls = list(),
	numsnp = list(),
	copynumbers = list(),
	conf = list(),
	tool.overlap.state = character(),
	overlap.coverage = list(),
	ref.tool = list(),
	n_genes = integer(),
	overlapping.genes = character()
	)
list_cols <- colnames(expected_final_tb)[sapply(expected_final_tb, function(x) is(x, 'list'))]

ensure_list_cols <- function(tb.or.gr){
	as_tibble(tb.or.gr) %>%
		rowwise() %>%
		mutate(across(one_of(list_cols), ~ list(.))) %>%
		makeGRangesFromDataFrame(keep.extra.columns = T)
}

overlap_tools <- function(tools, min.greater.region.overlap = 50) {
	
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
			reduce_ranges(tool = tool,
										individual_calls = individual_calls,
										numsnp = numsnp,
										copynumbers = copynumbers,
										conf = conf,
										widths = width,
										tool.overlap.state = ifelse(plyranges::n() > 1, 'post-overlap', 'no-overlap'),
										) %>%
			mutate( overlap.coverage =  round(widths / width * 100, 2),
							max.ov = sapply(overlap.coverage, max) ) %>%
			ensure_list_cols()
		ov[ov$tool.overlap.state == 'no-overlap',]$overlap.coverage <- NA
	} else {
		ov <- gr %>% 
			mutate(tool.overlap.state = 'no-overlap',
						 widths = NA,
						 overlap.coverage = NA,
						 max.ov = NA) %>%
			ensure_list_cols()
	}

	#Now filter based on (reciprocal) overlap to check if the overlap can be accepted
	# TODO: does the max(overlap.coverage) [i.e. not reciprocal] make sense?:
	# - reciprocal is hard to asses (we want to retain & check against final?)
	# - could only accept overlap if *all* tools have at least X% overlap with final (prevents ov chains)
	gr.not.changed <- ov %>% 
		filter(tool.overlap.state == 'no-overlap' | max.ov < min.greater.region.overlap) %>%
		plyranges::select(-widths, -max.ov) 
	
	gr.changed <- ov %>%
		filter(tool.overlap.state == 'post-overlap' & max.ov >= min.greater.region.overlap) %>%
		plyranges::select(-widths, -max.ov) 
		
	gr.pre.overlap <- find_overlaps(gr, gr.changed, suffix = c('', '.y')) %>% 
		filter(CNV.state == CNV.state.y) %>%
		plyranges::select(-contains('.y')) %>%
		mutate(tool.overlap.state = 'pre-overlap',
					 overlap.coverage = list(NA)) %>%
		ensure_list_cols()
		
	return(bind_ranges(gr.changed, gr.not.changed, gr.pre.overlap))
	
}


annotate_ref_overlap <- function(gr_in, gr_ref, min.cnv.coverage.by.ref = 0.8) {
	
	gr <- pair_overlaps(gr_in, gr_ref) 
	
	if (nrow(gr) > 0) {
		gr <- gr %>%
			as_tibble() %>% rowwise() %>%
			mutate(width.combined = max(granges.x.end, granges.y.end) - min(granges.x.start, granges.y.start) + 1,
					 width.overlap = min(granges.x.end, granges.y.end) - max(granges.x.start, granges.y.start) + 1,
					 coverage.by.ref = round(100 * width.overlap / granges.x.width, 2),
					 cov.ref.by.sample = round(100 * width.overlap / granges.y.width, 2),
					 ref.tool = list(tool.y),
					 ref.state = CNV.state.y
			) %>%
			dplyr::select(-contains('.y'), -contains('width')) %>%
			filter(coverage.by.ref >= min.cnv.coverage.by.ref & 
						 	CNV.state.x == ref.state) %>%
			group_by(granges.x.seqnames, granges.x.start, granges.x.end, tool.overlap.state.x, tool.x, CNV.state.x) %>%
			summarise(across(ends_with('.x'), ~ unique(.)),
								call.in.reference = TRUE,
								coverage.by.ref = list(coverage.by.ref),
								ref.tool = list(ref.tool),
			) %>%
			dplyr::rename_with(~ str_remove(., '.x$') %>% str_remove('^granges.x.'))

		# Need to check again, since checking for CNV state match can loose all matches, which breaks GRanges
		if (nrow(gr)>0){
		  gr <- makeGRangesFromDataFrame(gr, keep.extra.columns = T)
		} else {
		  gr <- GRanges()
		}


		# This keeps only (certain) overlaps, now
		# need to rebuild result if matching ref calls were found
		#1) Original regions not overlapping any merged call
		non.ovs <- filter_by_non_overlaps(gr_in, gr)
		# mutate fails on empty GRanges object
		if (length(non.ovs) > 0) { non.ovs <- non.ovs %>% mutate(call.in.reference = FALSE) }

		#2) Original regions overlapping that do overlap merged call (=in ref), but don't match the same call state.
		ovs.noStateMatch <- group_by_overlaps(gr_in, gr)
		if (length(ovs.noStateMatch) > 0) {
			ovs.noStateMatch <- ovs.noStateMatch %>%
				group_by(query, CNV.state.query) %>%
				mutate(any_match = any(CNV.state.query %in% CNV.state.subject)) %>%
				filter(!any_match) %>%
				as_tibble() %>%
				dplyr::select(-contains('.subject'), -any_match, -query) %>%
				dplyr::rename_with(~str_remove(., '.query'), contains('.query')) %>%
				mutate(call.in.reference = FALSE,
							 coverage.by.ref = list(NA),
							 ref.tool = list(NA)) %>%
				makeGRangesFromDataFrame(keep.extra.columns = T)
		}

		gr_out <- bind_ranges(
			# Calls with matching reference
			gr,
			# Calls not overlapping _any_ reference
			non.ovs,
			# Call overlapping a reference, but without any matching CNV.state
			ovs.noStateMatch
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

pennCNVfiles <- c(paste0(sample_id, '/', sample_id, '.penncnv-autosomes.tsv'),
									paste0(sample_id, '/', sample_id, '.penncnv-chrx.tsv'))
if (sex == 'm') { pennCNVfiles <- c(pennCNVfiles, paste0(sample_id, '/', sample_id, '.penncnv-chry.tsv')) }

results = list()

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

overlapped_tools_sample <- overlap_tools(tools, tool.overlap.min.perc)
 
if (!is.na(ref_id)) {
	
	fname <- file.path(datapath, ref_id, paste0(ref_id, '.combined-cnv-calls.tsv'))
	overlapped_tools_ref <- read_tsv(fname) %>%
		# These cols will interefere
		dplyr::select(-length, -call.in.reference, -coverage.by.ref, -ref.tool,
									-n_genes, -overlapping.genes) %>%
		mutate(across(one_of(list_cols), ~ str_split(., ';'))) %>%
		filter(tool.overlap.state != 'pre-overlap') %>%
		makeGRangesFromDataFrame(keep.extra.columns = T) 
	
	cnvs <- annotate_ref_overlap(overlapped_tools_sample, overlapped_tools_ref) %>%
		finalise_tb()
	
} else {
	cnvs <- overlapped_tools_sample %>% 
		finalise_tb()	
}
		
outname.tsv <- file.path(datapath, sample_id, paste0(sample_id, '.combined-cnv-calls.tsv'))
cnvs.tb <- cnvs %>%
	rowwise() %>%
	mutate(across(one_of(list_cols), ~paste(., collapse=';')))
write_tsv(cnvs.tb, outname.tsv)

#TODO add standard vcf output here instead of it's own script?
# -> can run another script to get different vcf, maybe?
