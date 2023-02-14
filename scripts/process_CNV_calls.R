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
args <- parse_args(parser, positional_arguments = 4, c("-p", "-c", "config.yaml", "sample_table_example.txt", "206210670080_R09C02", "data"))

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
tool.overlap.order <- config_val('tool.overlap.order', c('PennCNV', 'CBS'))
tool.overlap.min.perc <- config_val('tool.overlap.min.perc', 50)
# Compare to reference
min.cnv.coverage.by.ref <- config_val('min.cnv.coverage.by.ref', 80)

valid_name=config$wildcard_constraints$sample_id
valid_name=ifelse(is.null(valid_name), '[0-9]{12}_R[0-9]{2}C[0-9]{2}', valid_name)
if (!str_detect(sample_id, valid_name)) {stop('Sample id does not match supplied or default wildcard constraints!')}

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
			#TODO: does Y need adjustment as well?
			CNV.state = ifelse(Chr == 'chrX', NA, CNV.state),
			CNV.state = ifelse(Chr == 'chrX' & seg.median < -CBS.LRR.th.value + ifelse(sex == 'm', 1, -1) * CBS.LRR.th.value.Xadj, 'loss', CNV.state),
			CNV.state = ifelse(Chr == 'chrX' & seg.median > CBS.LRR.th.value + ifelse(sex == 'm', 1, -1) * CBS.LRR.th.value.Xadj, 'gain', CNV.state),
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

## Overlap tools
overlap_tools <- function(tool1.or.ov, tool2, min.greater.region.overlap = 50) {
	
	gr2 <- results[[tool2]] %>%
		mutate(tool.overlap.state = 'no-overlap')
	if (typeof(tool1.or.ov) == 'character') {
		gr1 <- results[[tool1.or.ov]] %>%
			mutate(tool.overlap.state = 'no-overlap',
						 overlap.coverage = NA,
						 overlapping.region = NA)
	} else {
		gr1 <- tool1.or.ov
	}
	
	tb.ovs <- gr1 %>%
		filter(tool.overlap.state != 'pre-overlap') %>%
		pair_overlaps(gr2) %>% 
		as_tibble() %>%
		rowwise() %>%
		mutate(
			#Keep size/region of call from first tool
			seqnames = granges.x.seqnames,
			start = granges.x.start,
			end = granges.x.end,
			sample_id = sample_id.x,
			#actual overlapping region
			#TODO: this will break shit if 3 intervals may have *NO* shared overlap -> need to find & remove somehow?
			start.ov = max(granges.x.start, granges.y.start, 
										 ifelse(!is.na(overlapping.region), as.integer(str_extract(overlapping.region, '(?<=:)[0-9]+(?=-)')), granges.x.start)),
			end.ov = min(granges.x.end, granges.y.end,
									 ifelse(!is.na(overlapping.region), as.integer(str_extract(overlapping.region, '(?<=-)[0-9]+$')), granges.x.end)),
			width.ov = end.ov -  start.ov + 1,
			#meta_cols
			#CNV.state = CNV.state.x,
			CNV.state = list(c(CNV.state.x, CNV.state.y)),
			tool = list(c(tool.x, tool.y)),
			#tool.calls = list(c(CNV.state.x, CNV.state.y)),
			#copynumbers = list(c(copynumbers.x, copynumbers.y)),
			numsnp = list(c(numsnp.x, numsnp.y)),
			overlapping.region = paste0(seqnames, ':', start.ov, '-', end.ov),
			individual_calls = list(c(individual_calls.x, individual_calls.y)),
			conf = list(c(conf.x, conf.y)),
			tool.overlap.state = 'post-overlap'
		) 
	
	# If coverage does not yet exist need to compute it for x (tool1)
	#  if tool1 was already an overlap the overlap.coverage should already exist
	tb.ovs[is.na(tb.ovs$overlap.coverage),]$overlap.coverage <- 
		with(tb.ovs[is.na(tb.ovs$overlap.coverage),], round(100 * width.ov / granges.x.width, 2)) %>%
		as.list()
	# Add coverage for y (tool2)
	tb.ovs <- tb.ovs %>% rowwise() %>%
		mutate(overlap.coverage = list(c(overlap.coverage, round(100 * width.ov / granges.y.width, 2))))		

	
	meta_cols = c('CNV.state', 'tool', 
								#'tool.calls', 
								'numsnp', 
								'overlap.coverage', 'overlapping.region', 'individual_calls','conf', 
								'tool.overlap.state')
	
	gr.ovs <- tb.ovs %>%
		dplyr::select(seqnames, start, end, sample_id, one_of(meta_cols)) %>%
		unique() %>%
		#Filter based on (reciprocal) overlap to check if the overlap is accepted?
		# TODO: does the max(overlap.coverage) [i.e. not quite reciprocal] make sense?
		# maybe use either reciprocal OR the other + match of CNV.state
		# --> then need special consideration of other matches with wrong CNV.state
		dplyr::filter(max(overlap.coverage) >= min.greater.region.overlap) %>%
		#TODO how to deal with i.e. LOH - loss/gain overlaps?
		# leave as is & extract within report?
		# > dont want to not overlap either, since the secondary (CBS) loss/gain call is more likely bullshit
		# could take first CNV.state or the one with better coverage / numnsp ?
		makeGRangesFromDataFrame(keep.extra.columns = T)
			
	gr.not.touched <- bind_ranges(gr1, gr2) %>% 
		filter(tool.overlap.state != 'pre-overlap') %>%
		filter_by_non_overlaps(gr.ovs)
	
	#TODO ? maybe somehow mark individual calls that don't match sate of higher order tools
	gr.overlapped <- bind_ranges(gr1, gr2) %>% 
		filter(tool.overlap.state == 'no-overlap') %>%
		filter_by_overlaps(gr.ovs) %>%
		mutate(tool.overlap.state = 'pre-overlap') %>%
		bind_ranges(filter(gr1, tool.overlap.state == 'pre-overlap'))
	
	bind_ranges(gr.ovs, gr.not.touched, gr.overlapped)
	
}


############
# Run code #
############

pennCNVfiles <- c(paste0(sample_id, '/', sample_id, '.penncnv-autosomes.', use.filter, '.tsv'),
									paste0(sample_id, '/', sample_id, '.penncnv-chrx.', use.filter, '.tsv'))
if (sex == 'm') { pennCNVfiles <- c(pennCNVfiles, paste0(sample_id, '/', sample_id, '.penncnv-chry.', use.filter, '.tsv')) }
if (!is.na(ref_id)) {
	pennCNVfiles_ref <- c(paste0(ref_id, '/', ref_id, '.penncnv-autosomes.', use.filter, '.tsv'),
												paste0(ref_id, '/', ref_id, '.penncnv-chrx.', use.filter, '.tsv'))
	if (sex == 'm') { pennCNVfiles_ref <- c(pennCNVfiles_ref, paste0(ref_id, '/', ref_id, '.penncnv-chry.', use.filter, '.tsv')) }
}


results = list()

if (args$options$penncnv) {
	results[['PennCNV']] <- file.path(datapath, pennCNVfiles) %>%
		lapply(read_PennCNV) %>%
		bind_rows()
	if (!is.na(ref_id)){
		results[['ref_PennCNV']] <- file.path(datapath, pennCNVfiles_ref) %>%
			lapply(read_PennCNV) %>%
			bind_rows()
	}
}
if (args$options$cbs) {
	results[['CBS']] <- file.path(datapath, sample_id, paste0(sample_id, '.CBS.', use.filter, '.tsv')) %>%
		lapply(read_CBS) %>%
		bind_rows()
	if (!is.na(ref_id)){
		results[['ref_CBS']] <- file.path(datapath, ref_id, paste0(ref_id, '.CBS.', use.filter, '.tsv')) %>%
			lapply(read_CBS) %>%
			bind_rows()
	}
}


if (merge.before.filter) {
	results <- lapply(results, merge_calls) %>%
		lapply(prefilter_calls)
} else {
	results <- lapply(results, prefilter_calls) %>%
		lapply(merge_calls)
}


# check which calls overlap between tools, calculate coverage & convert cols to lost, so the ycan store values from each tool post-overlap
# > for plotting original calls are still needed; these are retained (overlap.state == 'pre-overlap')


tools <- names(results) %>% str_subset('^ref_', negate = T)
if (!all(tool.overlap.order %in% tools)) {
	quit('"tool.overlap.order" contains tools that are not available')
} else if (!all(tools %in% tool.overlap.order)) {
	warning(paste0('Not all available CNV calling tools are specified in "tool.overlap.order".\n',
								 'These tools will *not* be used: ', paste0(tools[!tools %in% tool.overlap.order], collapse = ', ')))
}


if (length(tool.overlap.order) > 1) {
	overlapped_tools_sample <- tool.overlap.order[1]
	for (i in 2:length(tool.overlap.order)) {
		overlapped_tools_sample <- overlap_tools(overlapped_tools_sample, tool.overlap.order[i],
																	tool.overlap.min.perc)
	}
} else {
	message('Only 1 CNV calling tools used, no overlaps')
	#TODO:
	# probably need to convert the meta_cols into lists, to assure things work downstream
}


#TODO:
# maybe don't use overlapped ref calls, but the individual calls?
# & then do group_by
gr <- pair_overlaps(
							filter(overlapped_tools_sample, tool.overlap.state != 'pre-overlap'),
							bind_ranges(results[paste0('ref_', tool.overlap.order)])
							) %>% 
	as_tibble() %>% rowwise() %>%
	mutate(width.combined = max(granges.x.end, granges.y.end) - min(granges.x.start, granges.y.start) + 1,
				 width.overlap = min(granges.x.end, granges.y.end) - max(granges.x.start, granges.y.start) + 1,
				 coverage.by.ref = round(100 * width.overlap / granges.x.width, 2),
				 cov.ref.by.sample = round(100 * width.overlap / granges.y.width, 2),
				 ref.tool = list(tool.y),
				 ref.state = CNV.state.y[1],
				 sample.state = CNV.state.x[1]
				 ) %>%
	dplyr::select(-contains('.y'), -contains('width')) %>%
	filter(coverage.by.ref >= min.cnv.coverage.by.ref & 
				 	sample.state == ref.state) %>%
	group_by(granges.x.seqnames, granges.x.start, granges.x.end) %>%
	summarise(across(ends_with('.x'), ~ unique(.)),
						#tool.calls = unique(tool.calls),
						across(contains('overlap'), ~ unique(.)),
						call.in.reference = TRUE,
						coverage.by.ref = list(coverage.by.ref),
						ref.tool = list(ref.tool),
						) %>%
	dplyr::rename_with(~ str_remove(., '.x$') %>% str_remove('^granges.x.')) %>%
	makeGRangesFromDataFrame(keep.extra.columns = T)

cnvs <- bind_ranges(
	gr,
	filter_by_non_overlaps(
		filter(overlapped_tools_sample, tool.overlap.state != 'pre-overlap'),
		gr ),
	filter(overlapped_tools_sample, tool.overlap.state == 'pre-overlap'),
	) %>%
	as_tibble() %>%
	dplyr::rename(length = width) %>%
	dplyr::select(-strand)


#TODO: the list cols make saving a tsv cumbersome, though a tsv might be useful for export?
outname <- file.path(datapath, sample_id, paste0(sample_id, '.combined-cnv-calls.', use.filter, '.rds'))
saveRDS(cnvs, outname)
