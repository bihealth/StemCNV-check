suppressMessages(library(plyranges))
suppressMessages(library(GenomicRanges))
suppressMessages(library(tidyverse))
suppressMessages(library(optparse))
suppressMessages(library(yaml))

parser <- OptionParser(
	usage = "usage: %prog [options] data_path sample_id /path/to/config.yaml /path/to/sampletable.txt  "
)

parser <- add_option(parser, c("-p", "--penncnv"), default = F,
										 action="store_true", help="Use PennCNV calls")
parser <- add_option(parser, c("-c", "--cbs"), default = F,
										 action="store_true", help="Use CBS calls")

args <- parse_args(parser, positional_arguments = 4)
#args <- parse_args(parser, positional_arguments = 4, c("-p", "-c",  "data", "206764550040_R06C01", "config.yaml", "sample_table.txt"))

##################
# Variable setup #
##################

datapath <- args$args[1]
sample_id <- args$args[2]
config <- read_yaml(args$args[3])
sampletable <- read_tsv(args$args[4], col_types = 'cccccc') 

use.filter <- config$settings$filter$`use-filterset`
sex <- sampletable[sampletable$Sample_ID == sample_id, ]$Sex %>%
	tolower() %>% substr(1,1)

ref_name <- sampletable[sampletable$Sample_ID == sample_id, ]$ReferenceSample
ref_id <- NA

if (!is.na(ref_name)){
	ref_id <- sampletable[sampletable$SampleName == ref_name, ]$Sample_ID
	sex.ref <- sampletable[sampletable$Sample_ID == ref_id, ]$Sex %>%
		tolower() %>% substr(1,1)
	if(sex.ref != sex) {
		stop('Sex of sample and reference does not match!')
	} 
}


#TODO: might be easier to use a patched temp config file?
# Report only reads the actual config file, not the default one, so we need default values here as well
config_val <- function(value, default=NULL, section = 'postprocessing') {
	val=config$settings[[section]][[value]]
	if (is.null(val)) {
		return(default)
	} else {
		return(val)
	}
}

# CBS gain/loss
CBS.LRR.th.value <- config_val('LRR.th.value', 0.2, 'CBS')
CBS.LRR.th.value.Xadj <- config_val('LRR.th.value.Xadj', 0.2, 'CBS')
# Call merging & pre-filtering
min.snp					<- config_val('min.snp', 5)
min.length			<- config_val('min.length', 100)
min.snp.density <- config_val('min.snp.density', 10) #in snps per Mb
min.snp.or.length	<- config_val('min.snp.or.length', 10) # no.snps needed if length < min.length
min.length.or.snp	<- config_val('min.length.or.snp', 2000) #  length needed if no.snps < min.snp
merge.distance  <- config_val('merge.distance', 500) # max. distance for merging of small nearby calls
merge.before.filter  <- config_val('merge.before.filter', TRUE) # order of merging and filtering
# Tool overlapping
tool.order <- config_val('tool.order', c('PennCNV', 'CBS'))
tool.overlap.min.perc <- config_val('tool.overlap.min.perc', 50)
# Compare to reference
min.cnv.coverage.by.ref <- config_val('min.cnv.coverage.by.ref', 80)

valid_name=config$wildcard_constraints$sample_id
valid_name=ifelse(is.null(valid_name), '[0-9]{12}_R[0-9]{2}C[0-9]{2}', valid_name)
if (!str_detect(sample_id, valid_name)) {stop('Sample id does not match supplied or default wildcard constraints!')}

use_chr <- paste0('chr', c(1:22, 'X', 'Y'))

gtf_file <- config$static_data$genome_gtf_file
gr_genes  <- read_gff(gtf_file, col_names = c('source', 'type', 'gene_id', 'gene_type', 'gene_name')) %>%
	filter(type == 'gene' & seqnames %in% use_chr) 

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
					 sample_id = str_extract(sample_id, valid_name),
		)
}

## CBS
read_CBS <- function(filename) {
	tb <- read_tsv(filename, show_col_types = F) %>%
		dplyr::rename(Chr = chrom, start = loc.start, end = loc.end,
									numsnp = num.mark, sample_id = ID) %>%
		mutate(
			sample_id = str_remove(sample_id, '^X'),
			length = end - start, # TODO: open / half open / +- 1 ??
			Chr = paste0('chr', Chr),
			Chr = factor(Chr, levels = c(paste0('chr', 1:22), 'chrX', 'chrY')),
			snp.density = numsnp / length * 1e6,
			CNV.state = ifelse(seg.median < -CBS.LRR.th.value, 'loss', NA),
			CNV.state = ifelse(seg.median > CBS.LRR.th.value, 'gain', CNV.state),
			CNV.state = ifelse(Chr == 'chrX', NA, CNV.state),
			CNV.state = ifelse(Chr == 'chrX' & seg.median < -CBS.LRR.th.value + ifelse(sex == 'm', -1, 1) * CBS.LRR.th.value.Xadj, 'loss', CNV.state), 
			CNV.state = ifelse(Chr == 'chrX' & seg.median > CBS.LRR.th.value + ifelse(sex == 'm', -1, 1) * CBS.LRR.th.value.Xadj, 'gain', CNV.state), 
			tool = 'CBS',
			conf = NA,
			#TODO could try double / 3x / ... cutoff for 0/4 CN ?
			# > 0 doesn't make sense though (would mean no detection == filtered out)
			# > 4 maybe? very hard to tell & probably irrelevant
			copynumber = ifelse(is.na(CNV.state), 2, 3),
			copynumber = ifelse(CNV.state == 'loss', 1, copynumber),
		) %>%
		filter(!is.na(CNV.state) & !is.na(Chr))
	if (sex == 'f') {
		tb <- filter(tb, Chr != 'chrY')
	}
	tb
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

#TODO: move this to another script?
# Harmonize output
expected_final_tb = tibble(
	sample_id = character(),
	seqnames = factor(c(), levels = use_chr),
	start = integer(),
	end = integer(),
	length = integer(),
	#main.state = character(),
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
list_cols = colnames(expected_final_tb)[sapply(expected_final_tb, function(x) is(x, 'list'))]

ensure_list_cols <- function(tb.or.gr){ 
	
	as_tibble(tb.or.gr) %>%
		rowwise() %>%
		mutate(	across(one_of(list_cols), ~ list(.)),
						#This needs to happen inside rowwise
						#main.state = CNV.state[which.min(match(tool, tool.order))],
						) %>%
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
			dplyr::rename_with(~ str_remove(., '.x$') %>% str_remove('^granges.x.')) %>%
			makeGRangesFromDataFrame(keep.extra.columns = T)
		
		#Rebuild result if matching ref calls were found
		gr_out <- bind_ranges(
			# Calls with matching reference
			gr,
			# Calls not overlapping _any_ reference
			filter_by_non_overlaps(
				overlapped_tools_sample,
				gr) %>%
				mutate(call.in.reference = FALSE),
			# Call overlapping a reference, but without any matching CNV.state
			group_by_overlaps(overlapped_tools_sample, gr) %>% 
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

pennCNVfiles <- c(paste0(sample_id, '/', sample_id, '.penncnv-autosomes.', use.filter, '.tsv'),
									paste0(sample_id, '/', sample_id, '.penncnv-chrx.', use.filter, '.tsv'))
if (sex == 'm') { pennCNVfiles <- c(pennCNVfiles, paste0(sample_id, '/', sample_id, '.penncnv-chry.', use.filter, '.tsv')) }

results = list()

if (args$options$penncnv) {
	results[['PennCNV']] <- file.path(datapath, pennCNVfiles) %>%
		lapply(read_PennCNV) %>%
		bind_rows()
}
if (args$options$cbs) {
	results[['CBS']] <- file.path(datapath, sample_id, paste0(sample_id, '.CBS.', use.filter, '.tsv')) %>%
		lapply(read_CBS) %>%
		bind_rows()
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
	
	fname = file.path(datapath, ref_id, paste0(ref_id, '.combined-cnv-calls.', use.filter, '.rds'))
	overlapped_tools_ref <- readRDS(fname) %>%
		# These cols will interefere
		dplyr::select(-length, -call.in.reference, -coverage.by.ref, -ref.tool,
									-n_genes, -overlapping.genes) %>%
		filter(tool.overlap.state != 'pre-overlap') %>%
		makeGRangesFromDataFrame(keep.extra.columns = T) 
	
	cnvs <- annotate_ref_overlap(overlapped_tools_sample, overlapped_tools_ref) %>%
		finalise_tb()
	
} else {
	cnvs <- overlapped_tools_sample %>% 
		finalise_tb()	
}
		
#TODO: the list cols make saving a tsv cumbersome, though a tsv might be useful for export?
outname <- file.path(datapath, sample_id, paste0(sample_id, '.combined-cnv-calls.', use.filter, '.rds'))
saveRDS(cnvs, outname)
