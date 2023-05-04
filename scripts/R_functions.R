

filter_functions <- list(
	'basic'   = function(GT, GC) GT > 0.15,
	'highGT'   = function(GT, GC) GT > 0.8,
	'highGC'   = function(GT, GC) GT > 0.15 & GC > 0.8,
	'full' = function(GT, GC) GT > 0.8 & GC > 0.8,
	'highGT+GC' = function(GT, GC) GT > 0.8 & GC > 0.8
)


load_annotate_cnv <- function(fname){
	read_tsv(fname) %>%
		rowwise() %>%
		mutate(across(one_of(list_cols), ~ str_split(., ';'))) %>%
		mutate(
			#main.state = CNV.state[which.min(match(tool, tool.overlap.order))],
			#call.in.reference = ifelse(is.na(call.in.reference), F, call.in.reference),
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