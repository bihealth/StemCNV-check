#General
suppressMessages(require(GenomeInfoDb))
suppressMessages(require(tidyverse))
suppressMessages(require(readxl))
suppressMessages(require(plyranges))
`%!in%` <- Negate(`%in%`)

get_defined_labels <- function(config) {
    label_def_file <- file.path(config$snakedir, 'control_files', 'label_name_definitions.yaml')
    return(yaml.load_file(label_def_file))
}


fix_CHROM_format <- function(gr.or.tb, target_style) {
    stopifnot("GRanges" %in% class(gr.or.tb) | "data.frame" %in% class(gr.or.tb))
    if ("GRanges" %in% class(gr.or.tb)) {
        seqlevelsStyle(gr.or.tb) <- target_style
        return(sortSeqlevels(gr.or.tb))
    } else {
        orig_cols <- colnames(gr.or.tb)
        gr <- as_granges(gr.or.tb)
        seqlevelsStyle(gr) <- target_style
        tb <- sortSeqlevels(gr) %>%
            as_tibble() %>%
            select(all_of(orig_cols))
        return(tb)
    }
}


get_sex_chroms <- function(target_style) {
    sex_chroms <- genomeStyles('Homo_sapiens') %>%
        filter(sex)
    stopifnot(target_style %in% colnames(sex_chroms))
    sex_chroms %>% pull(target_style)
  #   # Legacy option: Determine chrom style and return sex chroms
  #   stopifnot("GRanges" %in% class(tb.or.gr) | "data.frame" %in% class(tb.or.gr))
  #   if ("GRanges" %in% class(tb.or.gr)) {
  #       chroms <- tb.or.gr
  #   } else if(is.factor(tb.or.gr$seqnames)) {
  #       chroms <- levels(tb.or.gr$seqnames)
  #   } else {
  #       chroms <- as.character(tb.or.gr$seqnames)
  #       stopifnot(length(chroms), "Can not determine chrom style from empty list!")
  #   }
  #   
  # sex_chroms %>%
  #       pull(chroms %>% seqlevelsStyle() %>% head(1))
}


read_sampletable <- function(filename, col_remove_regex = NA) {
    
    if (str_detect(filename, '\\.(txt|tsv)$')) {
        tb <- read_tsv(filename, comment = '#', show_col_types = F) 
    } else if(str_detect(filename, '\\.csv$')) {
        tb <- read_csv(filename, comment = '#', show_col_types = F) 
    } else if (str_detect(filename, '\\.xlsx$')) {
        tb <- read_excel(filename) %>%
            filter(!if_all(everything(), ~is.na(.))) %>%
            # read_excel doesn't support the 'comment' flag
            # need to recreate that behaviour while properly handling (NA) an empty first column
            filter(if_any(1, ~!str_detect(
                ifelse(is.na(.), '', as.character(.)),
                '^#')
            ))
    } else stop(paste('Unsupported file format:', filename))
    
    # Optional removal/editing of columns with a regex
    if (!is.na(col_remove_regex) & col_remove_regex != '') {
        tb <- rename_with(tb, ~str_remove(., col_remove_regex))
    }
    # Ensure all columns are characters
    mutate(tb, across(everything(), ~as.character(.)))
}


get_sample_info <- function(sample_id, value, config, sampletable = NA) {
    if ('data.frame' %in% class(sampletable)) { 
        sampletable <- sampletable 
    } else { 
        sampletable <- read_sampletable(config$sample_table, config$column_remove_regex)
    }
    # Handle NA input
    if (is.na(sample_id)) {
        return(NA_character_)
    }
    value_mapping <- c(
        'ref_id' = 'Reference_Sample',
        'sex' = 'Sex'
    )
    if (value %!in% c(colnames(sampletable), names(value_mapping))) {
        stop(paste('Unsupported sample info value:', value))
    }
    # remap legacy value names, then extract from sampletable
    if (value %in% names(value_mapping)) {
        value <- value_mapping[[value]]
    }

    out_val <- sampletable %>%
        filter(Sample_ID == sample_id) %>%
        pull(!!sym(value))
    # Coerce Sex to single character
    if (value == 'Sex') {
        out_val <- out_val %>% tolower() %>% substr(1, 1)
    }
    return(out_val)
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

load_gtf_data <- function(gtf_file, config, target_style='UCSC', include_hotspot_genes=c()) {
	exclude_regexes <- config$settings$CNV_processing$gene_overlap$exclude_gene_type_regex %>%
			paste(collapse = '|')
	gene_type_whitelist <- config$settings$CNV_processing$gene_overlap$include_only_these_gene_types
    
    whitelist_hotspots <- config$settings$CNV_processing$gene_overlap$whitelist_hotspot_genes
    if (whitelist_hotspots & length(include_hotspot_genes) == 0) {
        warning('No hotspot genes provided, ne genes will be whitelisted!')
    } else if (!whitelist_hotspots) {
        include_hotspot_genes <- c()   
    }

	gr_genes  <- read_gff(
        fix_rel_filepath(gtf_file, config),
        col_names = c('source', 'type', 'gene_id', 'gene_type', 'gene_name')
    ) %>%
		filter(type == 'gene') %>%
		mutate(gene_id = str_remove(gene_id, '\\..*')) %>%
        fix_CHROM_format(target_style)
    
	if (exclude_regexes != ''){
		gr_genes <- filter(gr_genes, !str_detect(gene_type, exclude_regexes) | gene_name %in% include_hotspot_genes)
	}
	if (is.character(gene_type_whitelist) & length(gene_type_whitelist) > 0){
		gr_genes <- filter(gr_genes, gene_type %in% gene_type_whitelist | gene_name %in% include_hotspot_genes)
	}
	gr_genes
}

## GenomeInfo Data

load_genomeInfo <- function(ginfo_file, config, target_style='UCSC') {
	# cols: chr	size	band_start	band_end	band_name	band_staining	centromer
	gr_info <- read_tsv(
            fix_rel_filepath(ginfo_file, config),
	        show_col_types = FALSE
        ) %>%
		filter(!is.na(band_start)) %>%
		as_granges(seqnames = chr, start = band_start, end = band_end) %>%
		mutate(section_name = paste0(str_remove(as.character(seqnames), 'chr'), band_name)) %>%
        fix_CHROM_format(target_style)

	gr_info
}


unsplit_merged_CNV_callers <- function(cnv_gr, defined_labels) {
    
    unmerged_calls <- cnv_gr %>%
        filter(CNV_caller != defined_labels$combined_cnvs) %>%
    
    merged_calls <- cnv_gr %>%
        filter(CNV_caller == defined_labels$combined_cnvs) %>%
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