
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
		mutate(across(any_of(list_cols), ~ str_split(., ';')))
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

get_chromosome_set <- function(use.config = NULL, get_prefix = F) {
	# Use given config if possible
	if(!is.null(use.config)) {
		use_chromosomes <- use.config$settings$chromosomes
	# else check outer scope for config object
	} else if(exists('config')) {
		use_chromosomes <- config$settings$chromosomes
	} else {
		warning('Setting Chromosomes to default wihtout using config!')
		use_chromosomes <- paste0('chr', c(1:22, 'X', 'Y'))
	}

	if (get_prefix) {
		chr_prefix <- ifelse(all(str_detect(use_chromosomes, '^chr')),
		                      'chr', '')
		return(chr_prefix)

	} else {
		return(use_chromosomes)
	}
}


# Harmonize output
expected_final_tb <- tibble(
	sample_id = character(),
	seqnames = factor(c(), levels = get_chromosome_set()),
	start = integer(),
	end = integer(),
	#Maybe make this a list_col ? granges can calculate width anyway
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


load_gtf_data <- function(config) {
	gtf_file <- config$static_data$genome_gtf_file
	#May need to ensure path is absolute, since rmd might have changed wd ?
	if (!file.exists(gtf_file)) {
		gtf_file <- normalizePath(gtf_file, mustWork = TRUE)
	}
	exclude_regexes <- config$settings$gene_overlap$exclude_gene_type_regex %>%
			paste(collapse = '|')
	gene_type_whitelist <- config$settings$gene_overlap$include_only_these_gene_types

	gr_genes  <- read_gff(gtf_file, col_names = c('source', 'type', 'gene_id', 'gene_type', 'gene_name')) %>%
		filter(type == 'gene' & !str_detect(gene_type, exclude_regexes))
	if (typeof(gene_type_whitelist) == 'list' & length(gene_type_whitelist) > 0){
		gr_genes <- filter(gr_genes, gene_type %in% gene_type_whitelist)
	}
	gr_genes
}

get_sample_info <- function(sample_id, value, sampletable){

	ref_id <- sampletable[sampletable$Sample_ID == sample_id, ]$Reference_Sample
	if (value == 'ref_id') return(ref_id)

	sex <- sampletable[sampletable$Sample_ID == sample_id, ]$Sex %>%
		tolower() %>% substr(1,1)
	if (value == 'sex') return(sex)

	if (!is.na(ref_id)){
		sex.ref <- sampletable[sampletable$Sample_ID == ref_id, ]$Sex %>%
  			tolower() %>% substr(1,1)
		if(sex.ref != sex) {
			stop('Sex of sample and reference does not match!')
		}
	} else {
		sex.ref <- NA
	}
	if (value == 'sex.ref') return(sex.ref)
	else stop(paste('Unsupported sample info value:', value))
}