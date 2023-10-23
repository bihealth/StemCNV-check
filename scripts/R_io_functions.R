#General
`%!in%` <- Negate(`%in%`)

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



# File input functions

## SNP data
read_raw <- function(filename) {
    read_tsv(filename, show_col_types = FALSE) %>%
                rename_with(~ str_remove(., '.*\\.')) %>%
        dplyr::select(-any_of(c('Index', 'Address', 'R', 'Theta')),
                                    -contains('Frac'), -contains('X'), -contains('Y')) %>%
        mutate(sample_id = basename(filename) %>% str_remove('\\.(processed|filtered)-data.*\\.tsv$'),
               Chr = ifelse(!str_detect(Chr, 'chr'), paste0('chr', Chr), Chr),
               )
}

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

## preprocessed

load_preprocessed_cnvs <- function(fname){
	read_tsv(fname) %>%
		rowwise() %>%
		mutate(across(any_of(list_cols), ~ str_split(., ';')))
}

## GTF data
get_static_path <- function(path, project_base=''){
	#Rmd might change cwd, so relative paths can break if not read/forwarded by snakemake
	if(file.exists(path)) return(path)
	else if (file.exists(file.path(project_base, path))) return(file.path(project_base, path))
	else if (file.exists(normalizePath(path))) return(normalizePath(path))
	else stop(paste('Could not find file path:', path))
}

load_gtf_data <- function(config) {
	gtf_file <- get_static_path(config$static_data$genome_gtf_file, config$basedir)
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


# Output

## Default table structure for CNVs
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

# Other report functions

get_SNP_clustering_IDs <- function(config_entry, sample_id, sampletable) {
    # '__[column]' entry: take all sample_ids with the same value in '[column]'
    same_value_cols <- str_subset(config_entry, '^__') %>% str_remove('^__')
    ids <- sapply(same_value_cols, function(x) {
        sample_value <- unlist(subset(sampletable, Sample_ID == sample_id)[, x], use.names = F)
        sampletable %>% filter(!!sym(x) == sample_value) %>% pull(Sample_ID)
    }) %>% unlist()
    # '_[column]' entry: take all sample_ids from '[column]'
    id_cols <- str_subset(config_entry, '^_[^_]') %>% str_remove('^_')
    ids <- c(ids, sapply(id_cols, function(x) {
        sampletable %>% filter(Sample_ID == sample_id) %>% pull(!!sym(x)) %>% str_split(',')
    }) %>% unlist() )
    # other entries: assume they are sample_ids & take the existing ones. No warning since this is discouraged
    single.ids <- str_subset(config_entry, '^[^_]')
    ids <- c(ids, single.ids[single.ids %in% sampletable$Sample_ID])
	ids
}