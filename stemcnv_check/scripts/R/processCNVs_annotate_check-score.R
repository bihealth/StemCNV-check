# Sanitize output & add gene overlap annotation
finalise_gr_to_tb <- function(gr, gr_genes) {

	gr$n_genes <- count_overlaps(gr, gr_genes)
	gr$overlapping_genes <- NA_character_
	ov_genes <- group_by_overlaps(gr, gr_genes) %>% reduce_ranges(genes = paste(gene_name, collapse = ','))
	gr[ov_genes$query,]$overlapping_genes <- ov_genes$genes

	tb <- as_tibble(gr) %>%
		rowwise() %>%
		mutate(across(any_of(list_cols), ~ list(.))) %>%
		bind_rows(expected_final_tb) %>%
		rowwise() %>%
		mutate(length = ifelse('lenght' %in% names(.), length, width)) %>%
		ungroup %>%
		dplyr::select(one_of(colnames(expected_final_tb)))

	return(tb)
}


annotate_cnv.check.score <- function(tb, high_impact_gr, highlight_gr, check_scores) {

	tb %>%
	  rowwise() %>%
	  mutate(
		  `Check-Score` =
			ifelse(CNV_type %in% c('gain', 'loss'),
				   1/3 * log(length) * log(length) - 15,
				   0.275 * log(length) * log(length) - 15
			) +
			# Base score for hitting HI / HL / ROI
			(check_scores$highimpact_base * !is.na(high_impact_hits) ) +
		    (check_scores$highlight_base * !is.na(highlight_hits) ) +
		    (check_scores$roi_hit_base * !is.na(ROI_hits) ) +
		    # Per gene/region score:
			# - scores per ROI (gene & others handeled the same)
			ifelse(!is.na(ROI_hits), (1 + str_count(ROI_hits, ',')) * check_scores$per_gene_roi, 0) +
			# - scores per non-gene HI / HL
			(high_impact_gr %>% as_tibble() %>%
				filter(mapping != 'gene_name' & hotspot %in% unlist(str_split(high_impact_hits, ','))) %>%
				mutate(check_score = ifelse(is.na(check_score), check_scores$per_gene_highimpact, check_score)) %>%
				pull(check_score) %>%	sum()) +
			(highlight_gr %>% as_tibble() %>%
				filter(mapping != 'gene_name' & hotspot %in% unlist(str_split(highlight_hits, ','))) %>%
				mutate(check_score = ifelse(is.na(check_score), check_scores$per_gene_highlight, check_score)) %>%
				pull(check_score) %>%	sum()) +
			# - score all genes that aren't also an ROI, use scores from tsv where available
			# dplyr will generate a bunch of warnings for calls without any genes
			suppressWarnings(bind_rows(
				high_impact_gr %>% as_tibble() %>%
					mutate(check_score = ifelse(is.na(check_score), check_scores$per_gene_highimpact, check_score)),
				highlight_gr %>% as_tibble() %>%
					mutate(check_score = ifelse(is.na(check_score), check_scores$per_gene_highlight, check_score)),
				str_split(overlapping_genes, ',') %>% unlist() %>%
					as_tibble() %>% dplyr::rename(hotspot = value) %>%
					mutate(mapping = 'gene_name', check_score = check_scores$per_gene_any)
			  	) %>%
				filter(hotspot %in% unlist(str_split(overlapping_genes, ',')) &
						   hotspot %!in% unlist(str_split(ROI_hits, ',')) &
						   mapping == 'gene_name') %>%
				group_by(hotspot) %>% summarise(check_score = max(check_score)) %>%
				pull(check_score) %>% sum())
	  ) %>%
	  ungroup() 
}

annotate_precision.estimates <- function(tb, size_categories, precision_estimates) {
	tb %>%
		rowwise() %>%
		mutate(
			size_category = case_when(
				length >= size_categories$extreme.loh & CNV_type %!in% c('gain', 'loss') ~ 'extreme',
				length >= size_categories$extreme.cnv & CNV_type %in% c('gain', 'loss') ~ 'extreme',
				length >= size_categories$very.large.loh & CNV_type %!in% c('gain', 'loss') ~ 'very_large',
				length >= size_categories$very.large.cnv & CNV_type %in% c('gain', 'loss') ~ 'very_large',
				length >= size_categories$large.loh & CNV_type %!in% c('gain', 'loss') ~ 'large',
				length >= size_categories$large.cnv & CNV_type %in% c('gain', 'loss') ~ 'large',
				length >= size_categories$medium.loh & CNV_type %!in% c('gain', 'loss') ~ 'medium',
				length >= size_categories$medium.cnv & CNV_type %in% c('gain', 'loss') ~ 'medium',
				TRUE ~ 'small'
		  	),
			Precision_Estimate =
			  ifelse(CNV_type %in% c('gain', 'loss'),
				precision_estimates[[
					ifelse(caller_merging_state == 'combined', 'multiple_Callers', unlist(CNV_caller))]][[
					size_category]] +
				  (precision_estimates$Call_has_Gap * (probe_coverage_gap)) +
				  (precision_estimates$HighSNPDensity * (high_probe_density)),
				NA_real_
			  )
		) %>%
		ungroup() %>%
	  dplyr::select(-size_category)
}