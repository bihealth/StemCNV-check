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
			seqnames = str_remove(!!sym(colname),':.*'),
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
	filter_regex <- paste0('^(',
		paste(str_replace(unlist(tb[, colname]), fixed('.'), '\\.'), collapse='|'),
		')')

	gr.tb <- gr_info %>%
		filter(str_detect(section_name, filter_regex)) %>%
		mutate(!!colname := str_extract(section_name, filter_regex)) %>%
		group_by(!!sym(colname)) %>%
		reduce_ranges() %>%
		as_tibble()

	not_matched <- tb[unlist(tb[colname]) %!in% unlist(gr.tb[colname]), colname]
	if (any(unlist(tb[colname]) %!in% unlist(gr.tb[colname]))) {
		stop(paste('The following band_names could not be identified in the reference data:', paste(not_matched, collapse = ', ')))
	}

	left_join(gr.tb, tb, by = colname) %>%
		as_granges()
}

# Hotspot list based annotation
parse_hotspot_table <- function(tb, gr_genes, gr_info) {

	if (any(tb$mapping %!in% c('gene_name', 'position', 'gband'))) {
		stop('Invalid mapping types in hotspot table')
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
        source = character(),
        comment = character()
    )

	if (nrow(sub_tb_name) > 0) {
		gr_name <- gr_genes %>%
			select(-source, -type, -gene_id, -gene_type) %>%
			filter(gene_name %in% sub_tb_name$hotspot) %>%
			as_tibble() %>%
			dplyr::rename(hotspot = gene_name) %>%
			left_join(sub_tb_name) %>%
			as_granges()
		#message('parsed gene names')
	} else {
		gr_name <- empty_gr
	}

	if (nrow(sub_tb_pos) > 0) {
		gr_pos <- tb_to_gr_by_position(sub_tb_pos, 'hotspot')
		#message('parsed gene positions')
	} else {
		gr_pos <- empty_gr
	}

	if (nrow(sub_tb_gband) > 0) {
		gr_gband <- tb_to_gr_by_gband(sub_tb_gband, gr_info, 'hotspot')
	} else {
		gr_gband <- empty_gr
	}
	bind_ranges(empty_gr, gr_name, gr_pos, gr_gband)
}

annotate_impact_lists <- function(gr, hotspot_gr, list_name) {
	message('Annotation calls with gene lists')

	# Make an extra col listing all overlapping directly defined genes
	gr@elementMetadata[[paste0(list_name, '_hits')]] <- NA_character_
	ov <- group_by_overlaps(gr, hotspot_gr)
	if (length(ov) > 0) {
		ov_hits <- ov %>%
			mutate(type_check = str_detect(CNV_type, str_replace(call_type, '^any$', '.*'))) %>%
			filter(type_check) %>%
			reduce_ranges(hits = paste(unique(sort(hotspot)),collapse = '|'))
		gr[ov_hits$query,]@elementMetadata[[paste0(list_name, '_hits')]] <- ov_hits$hits
	}

	return(gr)
}

annotate_roi <- function(gr, sample_id, sampletable, gr_genes, gr_info) {

	if (! 'Regions_of_Interest' %in% colnames(sampletable)) {
		return(gr %>% mutate(ROI_hits = NA_character_))
	}

	message('Annotating calls with ROI')

	roi <- sampletable[sampletable$Sample_ID == sample_id, ]$Regions_of_Interest

	if (is.na(roi) | is.null(roi) | roi == '') {
		return(gr %>% mutate(ROI_hits = NA_character_))
	}

	roi_gr <- roi %>% str_split(';') %>% unlist() %>%
		as_tibble() %>% dplyr::rename(roi = value) %>%
		mutate(hotspot = str_remove(roi, '^.*\\|'),
			   mapping = case_when(
				   str_detect(roi, '(chr)?[0-9XY]{1,2}:[0-9.,]+-[0-9.,]+') ~ 'position',
				   str_detect(roi, '[0-9XY]{1,2}(p|q)[0-9.]+') ~ 'gband',
				   TRUE ~ 'gene_name'
			   ),
			   roi_name = str_extract(roi, '^[^|]+\\|') %>% str_remove('\\|'),
			   roi_name = case_when(
				   mapping == 'gene_name' ~ hotspot,
				   is.na(roi_name) ~ paste0('ROI_', seq_along(roi_name)),
				   TRUE  ~ roi_name
			   ),
		) %>%
		parse_hotspot_table(gr_genes, gr_info)


	gr$ROI_hits <- NA_character_
	ov <- group_by_overlaps(gr, roi_gr)
	if (length(ov) > 0) {
		ov_hits <- ov %>%
			reduce_ranges(hits = paste(unique(roi_name), collapse = '|'))
		gr[ov_hits$query,]$ROI_hits <- ov_hits$hits
	}

	gr
}