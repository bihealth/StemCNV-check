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


fix_rel_filepath <- function(path, config){
	#Rmd might change cwd, so relative paths can break if not read/forwarded by snakemake
	if(file.exists(path)) return(path)
	else if (file.exists(file.path(config$basedir, path))) return(file.path(config$basedir, path))
	else if (file.exists(normalizePath(path))) return(normalizePath(path))
	else stop(paste('Could not find file path:', path))
}

load_gtf_data <- function(config, target_style='UCSC') {
	gtf_file <- fix_rel_filepath(config$static_data$genome_gtf_file, config)
	exclude_regexes <- config$settings$CNV_processing$gene_overlap$exclude_gene_type_regex %>%
			paste(collapse = '|')
	gene_type_whitelist <- config$settings$CNV_processing$gene_overlap$include_only_these_gene_types

	gr_genes  <- read_gff(gtf_file, col_names = c('source', 'type', 'gene_id', 'gene_type', 'gene_name')) %>%
		filter(type == 'gene') %>%
		mutate(gene_id = str_remove(gene_id, '\\..*')) %>%
        fix_CHROM_format(target_style)
    
	if (exclude_regexes != ''){
		gr_genes <- filter(gr_genes, !str_detect(gene_type, exclude_regexes))
	}
	if (is.character(gene_type_whitelist) & length(gene_type_whitelist) > 0){
		gr_genes <- filter(gr_genes, gene_type %in% gene_type_whitelist)
	}
	gr_genes
}

## GenomeInfo Data

load_genomeInfo <- function(config, target_style='UCSC') {
    
	# cols: chr	size	band_start	band_end	band_name	band_staining	centromer
	gr_info <- read_tsv(fix_rel_filepath(config$static_data$genomeInfo_file, config),
	                    show_col_types = FALSE) %>%
		filter(!is.na(band_start)) %>%
		as_granges(seqnames = chr, start = band_start, end = band_end) %>%
		mutate(section_name = paste0(str_remove(as.character(seqnames), 'chr'), band_name)) %>%
        fix_CHROM_format(target_style)

	gr_info
}

load_hotspot_table <- function(config, table = 'HighImpact') {
    
    if (table == 'HighImpact') {
        filename <- config$settings$CNV_processing$gene_overlap$high_impact_list 
    } else if (table == 'Highlight') {
        filename <- config$settings$CNV_processing$gene_overlap$highlight_list
    } else {
        stop('Invalid table name')
    }
    
    tb <- str_replace(filename, '__inbuilt__', config$snakedir) %>%
        fix_rel_filepath(config) %>%
        read_tsv(show_col_types = FALSE) 
    
    # str_glue with argument injection only works properly with named arguments that aren't numbers
    description_html_pattern <- str_replace_all(
        tb$description,
        '([:,] ?)(.+?)\\{([0-9]+)\\}(?=, ?|\\\\n|$)',
        '\\1<a href="{a\\3}">\\2</a>'
    )
    tb$description_htmllinks <- map2_chr(
        tb$description_doi, description_html_pattern, 
        \(doi, pattern) {
            args <- doi %>% str_split(', ?') %>% unlist() %>% set_names(paste0('a', 1:length(.)))
            rlang::inject(str_glue(pattern, !!!args))
        }
    )
    
    tb
}

unsplit_merged_CNV_callers <- function(cnv_gr) {
    
    unmerged_calls <- cnv_gr %>%
        filter(CNV_caller != 'StemCNV-check')
    
    merged_calls <- cnv_gr %>%
        filter(CNV_caller == 'StemCNV-check') %>%
        as_tibble() %>%
        separate_rows(initial_call_details, sep = '\\|') %>%
        mutate(
            CNV_caller = str_extract(initial_call_details, '^[^_]+'),
            start = str_extract(initial_call_details, '(?<=_)[0-9]+(?=-)') %>% as.integer(),
            end = str_extract(initial_call_details, '(?<=-)[0-9]+(?=_CN)') %>% as.integer(),
            CN = str_extract(initial_call_details, '(?<=CN)[0-9]+(?=_cov)') %>% as.integer(),
            # overlap_merged_call = str_extract(initial_call_details, '(?<=cov)[0-9]+\\.[0-9]+')
        ) %>%
        select(-width) %>%
        as_granges()
    
    bind_ranges(merged_calls, unmerged_calls)

}