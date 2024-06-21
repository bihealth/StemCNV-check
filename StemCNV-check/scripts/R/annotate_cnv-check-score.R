annotate_cnv.check.score <- function(tb, high_impact_gr, highlight_gr,
							 score_thresholds, impact_scores, precision_extimates) {

	tb %>%
	  rowwise() %>%
	  mutate(
		  size_category = case_when(
			length >= score_thresholds$extreme.loh & CNV_type %!in% c('gain', 'loss') ~ 'extreme',
			length >= score_thresholds$extreme.cnv & CNV_type %in% c('gain', 'loss') ~ 'extreme',
			length >= score_thresholds$very.large.loh & CNV_type %!in% c('gain', 'loss') ~ 'very_large',
			length >= score_thresholds$very.large.cnv & CNV_type %in% c('gain', 'loss') ~ 'very_large',
			length >= score_thresholds$large.loh & CNV_type %!in% c('gain', 'loss') ~ 'large',
			length >= score_thresholds$large.cnv & CNV_type %in% c('gain', 'loss') ~ 'large',
			length >= score_thresholds$medium.loh & CNV_type %!in% c('gain', 'loss') ~ 'medium',
			length >= score_thresholds$medium.cnv & CNV_type %in% c('gain', 'loss') ~ 'medium',
			TRUE ~ 'small'
		  ),
		  Impact_Score =
			ifelse(CNV_type %in% c('gain', 'loss'),
				   1/3 * log(length) * log(length) - 15,
				   0.275 * log(length) * log(length) - 15
			) +
			# Base score for hitting HI / HL / ROI
			(impact_scores$highimpact_base * !is.na(high_impact_hits) ) +
		    (impact_scores$highlight_base * !is.na(highlight_hits) ) +
		    (impact_scores$roi_hit_base * !is.na(ROI_hits) ) +
		    # Per gene/region score:
			# - scores per ROI (gene & others handeled the same)
			ifelse(!is.na(ROI_hits), (1 + str_count(ROI_hits, ',')) * impact_scores$pere_gene_roi, 0) +
			# - scores per non-gene HI / HL
			(high_impact_gr %>% as_tibble() %>%
				filter(mapping != 'gene_name' & hotspot %in% unlist(str_split(high_impact_hits, ','))) %>%
				mutate(impact_score = ifelse(is.na(impact_score), impact_scores$per_gene_highimpact, impact_score)) %>%
				pull(impact_score) %>%	sum()) +
			(highlight_gr %>% as_tibble() %>%
				filter(mapping != 'gene_name' & hotspot %in% unlist(str_split(highlight_hits, ','))) %>%
				mutate(impact_score = ifelse(is.na(impact_score), impact_scores$per_gene_highlight, impact_score)) %>%
				pull(impact_score) %>%	sum()) +
			# - score all genes that aren't also an ROI, use scores from tsv where available
			# dplyr will generate a bunch of warnings for calls without any genes
			suppressWarnings(bind_rows(
				high_impact_gr %>% as_tibble() %>%
					mutate(impact_score = ifelse(is.na(impact_score), impact_scores$per_gene_highimpact, impact_score)),
				highlight_gr %>% as_tibble() %>%
					mutate(impact_score = ifelse(is.na(impact_score), impact_scores$per_gene_highlight, impact_score)),
				str_split(overlapping_genes, ',') %>% unlist() %>%
					as_tibble() %>% rename(hotspot = value) %>%
					mutate(mapping = 'gene_name', impact_score = impact_scores$per_gene_any)
			  	) %>%
				filter(hotspot %in% unlist(str_split(overlapping_genes, ',')) &
						   hotspot %!in% unlist(str_split(ROI_hits, ',')) &
						   mapping == 'gene_name') %>%
				group_by(hotspot) %>% summarise(impact_score = max(impact_score)) %>%
				pull(impact_score) %>% sum())
		  ,
		  Precision_Estimate =
			  ifelse(CNV_type %in% c('gain', 'loss'),
				precision_extimates[[
					ifelse(caller_merging_state == 'combined', 'multiple_Callers', unlist(CNV_caller))]][[
					size_category]] +
				  (precision_extimates$Call_has_Gap * (probe_coverage_gap)) +
				  (precision_extimates$HighSNPDensity * (high_probe_density)),
				NA_real_
			  )
	  ) %>%
	  ungroup() %>%
	  dplyr::select(-size_category)

}