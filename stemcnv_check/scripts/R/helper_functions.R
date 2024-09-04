#General
suppressMessages(require(GenomeInfoDb))
suppressMessages(require(tidyverse))
suppressMessages(require(plyranges))
`%!in%` <- Negate(`%in%`)

fix_CHROM_format <- function(gr, target_style) {
    seqlevelsStyle(gr) <- target_style
    return(sortSeqlevels(gr))
}

get_sex_chroms <- function(tb.or.gr) {
    if ("GRanges" %in% class(tb.or.gr)) {
        chroms <- tb.or.gr
    } else {
        chroms <- as.character(tb.or.gr$seqnames)
    }
    
    genomeStyles('Homo_sapiens') %>%
        filter(sex) %>%
        pull(chroms %>% seqlevelsStyle() %>% head(1))
}

read_sampletable <- function(filename) {
    read_tsv(filename, col_types = 'cccccc', comment = '#')
}


get_sample_info <- function(sample_id, value, sampletable) {
    if ('data.frame' %in% class(sampletable)) { sampletable <- sampletable }
    else if ('character' %in% class(sampletable)) { sampletable <- read_sampletable(sampletable) }
    else { stop('`sampletable` needs to be a table or file path') }

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

get_target_chrom_style <- function(config, snp_vcf_gr) {
    target_style <- config$settings$vcf_output$chrom_style
    if (target_style == 'keep-original') {
        target_style <- seqlevelsStyle(snp_vcf_gr) %>% head(1) 
    }
    return(target_style)
}

## preprocessed

load_preprocessed_cnvs <- function(fname){
	read_tsv(fname) %>%
		rowwise() %>%
		mutate(across(any_of(get_list_cols()), ~ str_split(., ';')))
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
	exclude_regexes <- config$settings$CNV_processing$gene_overlap$exclude_gene_type_regex %>%
			paste(collapse = '|')
	gene_type_whitelist <- config$settings$CNV_processing$gene_overlap$include_only_these_gene_types

	gr_genes  <- read_gff(gtf_file, col_names = c('source', 'type', 'gene_id', 'gene_type', 'gene_name')) %>%
		filter(type == 'gene') %>%
		mutate(gene_id = str_remove(gene_id, '\\..*'))
	if (exclude_regexes != ''){
		gr_genes <- filter(gr_genes, !str_detect(gene_type, exclude_regexes))
	}
	if (is.character(gene_type_whitelist) & length(gene_type_whitelist) > 0){
		gr_genes <- filter(gr_genes, gene_type %in% gene_type_whitelist)
	}
	gr_genes
}

## GenomeInfo Data

load_genomeInfo <- function(config) {

	# cols: chr	size	band_start	band_end	band_name	band_staining	centromer
	gr_info <- read_tsv(get_static_path(config$static_data$genomeInfo_file, config$basedir),
	                    show_col_types = FALSE) %>%
		filter(!is.na(band_start)) %>%
		as_granges(seqnames = chr, start = band_start, end = band_end) %>%
		mutate(section_name = paste0(str_remove(as.character(seqnames), 'chr'), band_name))

	gr_info
}

# Output

## Default table structure for CNVs
get_expected_final_tb <- function(chrom_style='UCSC') {

    tibble(
        sample_id = character(),
        seqnames = factor(c(), levels = genomeStyles('Homo_sapiens')[[chrom_style]]),
        start = integer(),
        end = integer(),
        width = integer(),
        CNV_type = character(),
        ID = character(),
        Check_Score = double(),
        reference_overlap = logical(), # not needed for vcf output
        CNV_caller = list(),
        # n_premerged_calls = list(),
        n_probes = integer(),
        n_uniq_probes = integer(),
        probe_density_Mb = double(),
        CN = integer(),
        LRR = double(),
        Precision_Estimate = double(),
        caller_merging_state = character(),
        overlap_merged_call = double(),
        caller_merging_coverage = character(),
        reference_caller = character(), # deprecate this
        reference_coverage = double(),
        high_impact_hits = character(),
        highlight_hits = character(),
        ROI_hits = character(),
        percent_gap_coverage = double(),
        probe_coverage_gap = logical(),
        high_probe_density = logical(),
        n_genes = integer(),
        overlapping_genes = character()
    )
}

get_list_cols <- function() {
    colnames(get_expected_final_tb())[sapply(get_expected_final_tb(), function(x) is(x, 'list'))]
}


ensure_list_cols <- function(tb.or.gr){
	as_tibble(tb.or.gr) %>%
		rowwise() %>%
		mutate(across(any_of(get_list_cols()), ~ list(.))) %>%
		as_granges()
}

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

