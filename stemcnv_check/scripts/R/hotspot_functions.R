suppressMessages(require(tidyverse))
suppressMessages(require(plyranges))

tb_to_gr_by_position <- function(tb, colname = 'position') {

	if(nrow(tb) == 0){
		tb %>%
			mutate(seqnames = character(),
			       start = integer(),
					end = integer()) %>%
			as_granges() %>%
			return()
	}

	no_match <- tb %>% filter(!str_detect(!!sym(colname), '(chr)?[0-9XY]{1,2}:[0-9.,]+-[0-9.,]+'))
	if (nrow(no_match) > 0) {
		stop(paste('The following positions can not be converted:', paste(no_match[, colname], collapse = ', ')))
	}

	gr <- tb %>%
		filter(str_detect(!!sym(colname), '(chr)?[0-9XY]{1,2}:[0-9.,]+-[0-9.,]+')) %>%
		mutate(
            # need to return consistent seqnames
			seqnames = str_remove(!!sym(colname),':.*') %>% str_remove('^chr'),
			start = str_extract(!!sym(colname), '(?<=:)[0-9]+') %>% as.numeric(),
			end = str_extract(!!sym(colname), '(?<=-)[0-9]+$') %>% as.numeric()
		) %>% as_granges()

	return(gr)
}

tb_to_gr_by_gband <- function(tb, gr_info, colname = 'band_name') {

	if(nrow(tb) == 0){
		return(tb %>%
			mutate(seqnames = character(),
			       start = integer(),
					end = integer()) %>%
			as_granges()
			)
	}

	no_match <- tb %>% filter(!str_detect(!!sym(colname), '[0-9XY]{1,2}(p|q)[0-9.]+')) %>%
		pull(!!sym(colname))
	if (length(no_match) > 0) {
		stop(paste('The following band_names can not be converted:', paste(no_match, collapse = ', ')))
	}

	tb <- tb %>% filter(str_detect(!!sym(colname), '[0-9XY]{1,2}(p|q)[0-9.]+'))

    # get a _set_ of matching bands (=gr rows) for each tb row
    # then do the str_extract and group_by on each set
    gr.tb <- paste0(
        '^', str_replace(unlist(tb[, colname]), fixed('.'), '\\.')
    ) %>%        
        lapply(\(gband_regex) {
            gr_info %>%
                filter(str_detect(section_name, gband_regex)) %>%
                mutate(!!colname := str_extract(section_name, gband_regex)) %>%
                group_by(!!sym(colname)) %>%
                reduce_ranges() %>%
                as_tibble()
        }) %>%
        bind_rows() %>%
        # discard multiples in case of same gband for i.e. loss & gain
        unique()

	not_matched <- tb[unlist(tb[colname]) %!in% unlist(gr.tb[colname]), colname]
	if (any(unlist(tb[colname]) %!in% unlist(gr.tb[colname]))) {
		stop(paste('The following band_names could not be identified in the reference data:', paste(not_matched, collapse = ', ')))
	}

	left_join(gr.tb, tb, by = colname) %>%
		as_granges()
}

# Hotspot list based annotation
parse_hotspot_table <- function(tb, gr_genes, gr_info, target_chrom_style=NA) {

	if (any(tb$mapping %!in% c('gene_name', 'position', 'gband'))) {
		stop('Invalid mapping types in hotspot table')
	}
    
    if (is.na(target_chrom_style) | target_chrom_style == 'keep-original') {
        target_chrom_style <- seqlevelsStyle(gr_genes)
    }

	sub_tb_name <- tb %>%
		filter(mapping == 'gene_name')
	sub_tb_pos <- tb %>%
		filter(mapping == 'position')
	sub_tb_gband <- tb %>%
		filter(mapping == 'gband')

	empty_gr <- GRanges(
        list_name = character(),
        hotspot = character(),
        mapping = character(),
        call_type = character(),
        check_score = numeric(),
        description = character(),
        description_doi = character()
    )

	if (nrow(sub_tb_name) > 0) {
		gr_name <- gr_genes %>%
			select(-source, -type, -gene_id, -gene_type) %>%
			filter(gene_name %in% sub_tb_name$hotspot) %>%
			as_tibble() %>%
			dplyr::rename(hotspot = gene_name) %>%
            mutate(call_type = 'any;loss;gain;LOH') %>%
            separate_rows(call_type, sep = ';') %>%
            # Specifically merge each call_type
			inner_join(sub_tb_name, by = c('hotspot', 'call_type')) %>%
            # Remove duplicates (can result from non-unique gene_names)
            slice_max(check_score, by = c(hotspot, call_type), with_ties = FALSE) %>%            
			as_granges() %>%
            fix_CHROM_format(target_chrom_style)
		#message('parsed gene names')
        unmatched_genes <- setdiff(sub_tb_name$hotspot, gr_name$hotspot)
        if (length(unmatched_genes) > 0) {
            message('The following gene names could not be identified in the gtf file (they are likely alternave names):', paste(unmatched_genes, collapse = ', '))
        }
        
	} else {
		gr_name <- empty_gr
	}

	if (nrow(sub_tb_pos) > 0) {
		gr_pos <- tb_to_gr_by_position(sub_tb_pos, 'hotspot') %>%
            fix_CHROM_format(target_chrom_style)
	} else {
		gr_pos <- empty_gr
	}

	if (nrow(sub_tb_gband) > 0) {
		gr_gband <- tb_to_gr_by_gband(sub_tb_gband, gr_info, 'hotspot') %>%
            fix_CHROM_format(target_chrom_style)
	} else {
		gr_gband <- empty_gr
	}
	bind_ranges(empty_gr, gr_name, gr_pos, gr_gband)
}


get_dosage_sensivity_tb <- function(score_settings) {
    
    dosage_data <- 'https://zenodo.org/records/6347673/files/Collins_rCNV_2022.dosage_sensitivity_scores.tsv.gz'
    doi <- '10.1016/j.cell.2022.06.036'
    
    tb <- read_tsv(dosage_data, show_col_types = F) %>%
        rename_with(
            ~str_replace(., '#gene', 'hotspot') %>%
                str_replace('pHaplo', 'loss') %>%
                str_replace('pTriplo', 'gain')        
        ) %>%
        pivot_longer(cols = -hotspot, names_to = 'call_type', values_to = 'dosage_score') %>%
        mutate(
            list_name = 'Dosage-sensivity',
            mapping = 'gene_name',
            check_score = case_when(
                call_type == 'loss' & dosage_score >= score_settings$pHaplo_threshold ~ score_settings$dosage_sensitive_gene,
                call_type == 'gain' & dosage_score >= score_settings$pTriplo_threshold ~ score_settings$dosage_sensitive_gene,
                TRUE ~ NA_integer_
            ),
            description = paste0(
                'Gene with predicted dosage sensitivity (',
                ifelse(call_type == 'loss', 'haploinsufficiency', 'triplosensitivity'), ')\\n',
                'Source: Collins et al. 2022 {1}.\\n',
                ifelse(call_type == 'loss', 'pHaplo', 'pTriplo'),
                ' score: ', round(dosage_score, 3)                
            ),            
            description_doi = doi,
        ) %>%
        filter(!is.na(check_score)) %>%
        select(-dosage_score)
    
    description_html_pattern <- str_replace_all(
        tb$description,
        '([:,] ?)([^:,]+?)\\{([0-9]+)\\}(?=, ?|\\\\\\\\n|\\\\n|$)',
        '\\1<a href="{a\\3}" target="_blank" rel="noopener noreferrer">\\2</a>'
    ) %>%
        str_replace_all('\\\\n', '&#013;') %>%
        str_replace_all('\\n', '&#013;')
    
    tb$description_htmllinks <- map2_chr(
        tb$description_doi, description_html_pattern, 
        \(doi, pattern) {
            args <- doi %>% 
                str_split(', ?') %>% 
                # Make the doi text into an actual link
                sapply(\(x) paste0('https://doi.org/', x)) %>%
                unlist() %>%
                set_names(paste0('a', 1:length(.)))
            rlang::inject(str_glue(pattern, !!!args))
        }
    )
    
    tb
}


load_hotspot_table <- function(config, table = 'stemcell_hotspot') {
    
    if (table == 'stemcell_hotspot') {
        filename <- config$settings$CNV_processing$gene_overlap$stemcell_hotspot_list 
    } else if (table == 'cancer_gene') {
        filename <- config$settings$CNV_processing$gene_overlap$cancer_gene_list
    } else if (table == 'snv_hotspot') {
        filename <- config$settings$SNV_analysis$snv_hotspot_table        
    } else {
        stop('Invalid table name')
    }
    
    tb <- str_replace(filename, '__inbuilt__', config$snakedir) %>%
        fix_rel_filepath(config) %>%
        read_tsv(show_col_types = FALSE) 
    
    # str_glue with argument injection only works properly with named arguments that aren't numbers
    description_html_pattern <- str_replace_all(
        tb$description,
        '([:,] ?)([^:,]+?)\\{([0-9]+)\\}(?=, ?|\\\\\\\\n|\\\\n|$)',
        '\\1<a href="{a\\3}" target="_blank" rel="noopener noreferrer">\\2</a>'
    ) %>%
        str_replace_all('\\\\n', '&#013;') %>%
        str_replace_all('\\n', '&#013;')
    tb$description_htmllinks <- map2_chr(
        tb$description_doi, description_html_pattern, 
        \(doi, pattern) {
            args <- doi %>% 
                str_split(', ?') %>% 
                # Make the doi text into an actual link
                sapply(\(x) paste0('https://doi.org/', x)) %>%
                unlist() %>%
                set_names(paste0('a', 1:length(.)))
            rlang::inject(str_glue(pattern, !!!args))
        }
    )
    
    if (table == 'snv_hotspot') {
        tb <- tb %>%
            mutate(
                gene_name = hotspot %>% str_remove('::.*'),
                HGVS.p = paste0(
                    'p.',
                    hotspot %>% str_remove('^.*::')
                )
            )
    }
    
    tb
}


get_roi_gr <- function(sample_id, sampletable, config, gr_genes, gr_info, target_chrom_style) {
    
    empty_res <- GRanges(
        list_name = character(),
        hotspot = character(),
        mapping = character(),
        call_type = character(),
        check_score = numeric(),
        description = character(),
        description_doi = character(),
        description_htmllinks = character()
    )
    if ('Regions_of_Interest' %!in% colnames(sampletable)) {
        return(empty_res)
    }
    
    roi_regions <- sampletable %>% 
        filter(Sample_ID == sample_id) %>% 
        pull(Regions_of_Interest) %>% 
        str_split(';') %>% unlist()
    
    if (all(roi_regions == '' | is.na(roi_regions))) {
        return(empty_res)
	} 
        
    tibble(
        list_name = 'ROI',
        description = str_extract(roi_regions, '^[^|]+\\|') %>% str_remove('\\|'),
        hotspot = str_remove(roi_regions, '^[^|]+\\|'),
        mapping = case_when(
            str_detect(hotspot, '^(chr)?[0-9XY]{1,2}[pq][0-9.]+') ~ 'gband',
            str_detect(hotspot, '^(chr)?[0-9XY]{1,2}:[0-9]+-[0-9]+') ~ 'position',
            TRUE ~ 'gene_name'
        ),
        call_type = 'any',
        check_score = config$settings$CNV_processing$Check_score_values$any_roi_hit,
        description_doi = NA_character_,
    ) %>%
        mutate(
            description = paste0(
                'ROI_', seq_along(description), ': ',
                hotspot,
                ifelse(
                    !is.na(description), 
                    paste(';', description),
                    ''
                )
            ),
            order = row_number(),
            description_htmllinks = description
        ) %>%
		parse_hotspot_table(gr_genes, gr_info) %>%
        fix_CHROM_format(target_chrom_style) %>%
        arrange(order) %>% select(-order)
}
