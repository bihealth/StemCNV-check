
`%!in%` <- Negate(`%in%`)

## PennCNV
read_PennCNV <- function(filename) {
	read.table(filename, sep='', header = F, fill=T,
						 col.names = c('Position', 'numsnp', 'length', 'hmm.state', 'input', 'startsnp', 'endsnp', 'tool_confidence')) %>%
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
			   ID = paste(tool, CNV.state, Chr, start, end, sep='_'),
			   # basename can't handly empty input ...
			   sample_id = str_remove(sample_id, '.*/') %>% str_remove('\\.filtered-data-.*\\.tsv$'),
		)
}

## CBS
read_CBS <- function(filename) {
	read_tsv(filename, show_col_types = F)
}


## GADA / MAD

read_GADA <- function(filename) {
	read_tsv(filename, show_col_types = F)
}


load_preprocessed_cnvs <- function(fname){
	read_tsv(fname) %>%
		rowwise() %>%
		mutate(across(one_of(list_cols), ~ str_split(., ';')))
}

annotate_cnvs <- function(tb.or.gr, config) {

	#Reportable thresholds
	reportable.loh <- config$settings$report$thresholds$reportable.loh %>% as.numeric()
	reportable.cnv <- config$settings$report$thresholds$reportable.cnv %>% as.numeric()
	fail.loh <- config$settings$report$thresholds$fail.loh %>% as.numeric()
	fail.cnv <- config$settings$report$thresholds$fail.cnv %>% as.numeric()

	tb.or.gr %>%
			mutate(
			reportable = case_when(
				#tool.overlap.state == 'pre-overlap'                                ~ NA_character_,
				CNV.state %!in% c('gain', 'loss') & length >= reportable.loh & call.in.reference ~ 'yes, in ref.',
				CNV.state %in% c('gain', 'loss') & length >= reportable.cnv & call.in.reference ~ 'yes, in ref.',
				CNV.state %!in% c('gain', 'loss') & length >= fail.loh                           ~ 'critical',
				CNV.state %in% c('gain', 'loss') & length >= fail.cnv                           ~ 'critical',
				CNV.state %!in% c('gain', 'loss') & length >= reportable.loh                     ~ 'yes',
				CNV.state %in% c('gain', 'loss') & length >= reportable.cnv                     ~ 'yes',
				call.in.reference                                                 ~ 'no, in ref.',
				TRUE                                                              ~ 'no') %>%
				factor(levels = c('critical', 'yes', 'yes, in ref.', 'no', 'no, in ref.'))
		)
}


use_chr <- paste0('chr', c(1:22, 'X', 'Y'))

#TODO: move this + input function defs to another script?
# Harmonize output
expected_final_tb <- tibble(
	sample_id = character(),
	seqnames = factor(c(), levels = use_chr),
	start = integer(),
	end = integer(),
	#Maybe make this a lis_col ? granges can calculate width anyway
	length = integer(),
	CNV.state = character(),
	ID = character(),
	call.in.reference = logical(),
	coverage.by.ref = list(),
	tool = list(),
	merged_tool_calls = list(),
	numsnp = list(),
	copynumber = list(),
	tool_confidence = list(),
	tool.overlap.state = character(),
	coverage.overlap = list(),
	tool.coverage.overlap = character(),
	ref.tool = list(),
	n_genes = integer(),
	overlapping.genes = character()
	)
list_cols <- colnames(expected_final_tb)[sapply(expected_final_tb, function(x) is(x, 'list'))]

ensure_list_cols <- function(tb.or.gr){
	as_tibble(tb.or.gr) %>%
		rowwise() %>%
		mutate(across(any_of(list_cols), ~ list(.))) %>%
		as_granges()
}


# Sanitize output & add gene overlap annotation
finalise_cnv_tb <- function(gr) {

	#TODO this section should be moved
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
		dplyr::select(any_of(colnames(expected_final_tb)))

	return(tb)
}