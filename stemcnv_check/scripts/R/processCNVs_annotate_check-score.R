suppressMessages(require(tidyverse))
suppressMessages(require(plyranges))

# Sanitize output & add gene overlap annotation
annotate_gene_overlaps <- function(gr, gr_genes) {

	gr$n_genes <- count_overlaps(gr, gr_genes)
	gr$overlapping_genes <- NA_character_
	ov_genes <- group_by_overlaps(gr, gr_genes) %>% reduce_ranges(genes = paste(gene_name, collapse = '|'))
	gr[ov_genes$query,]$overlapping_genes <- ov_genes$genes

	return(gr)
}


annotate_cnv.check.score <- function(tb, stemcell_hotspots_gr, dosage_sensitive_gene_gr, cancer_genes_gr, check_scores) {

	tb %>%
	  rowwise() %>%
	  mutate(
		  Check_Score =
			ifelse(
                #Note: adapt this if CNV is used for CN >= 4
                CNV_type %in% c('gain', 'loss'),
				1/3 * log(width) * log(width) - 15,
				0.275 * log(width) * log(width) - 15
			) * 
            # Size modifier for large copynumbers
            ifelse(
                CN < 1 | CN > 3,
                check_scores$large_CN_size_modifier,
                1
            ) +
			# Base score for any ROI hit
			ifelse(!is.na(ROI_hits), check_scores$any_roi_hit, 0) +
            # Scores for each hotspot hit/dosage gene/cancer gene (cumulative, even if its the same gene)
			(stemcell_hotspots_gr %>% as_tibble() %>%
				filter(hotspot %in% unlist(str_split(stemcell_hotspot, '\\|'))) %>%
				pull(check_score) %>%	sum()) + 
            # dosage_sensitive & cancer lists should only have genes, but are included anyway
            (dosage_sensitive_gene_gr %>% as_tibble() %>%
				filter(hotspot %in% unlist(str_split(dosage_sensitive_gene, '\\|'))) %>%
				pull(check_score) %>%	sum()) +    
			(cancer_genes_gr %>% as_tibble() %>%
				filter(hotspot %in% unlist(str_split(cancer_gene, '\\|'))) %>%
				pull(check_score) %>%	sum()) +
			# Score all remaining (not hotpsot/dosage/cancer and not ROI)
			# dplyr will generate a bunch of warnings for calls without any genes
            str_split(overlapping_genes, '\\|') %>% unlist() %>%
                as_tibble() %>% dplyr::rename(gene_name = value) %>%
                filter(
                    gene_name %!in% c(
                        unlist(str_split(stemcell_hotspot, '\\|')),
                        unlist(str_split(dosage_sensitive_gene, '\\|')),
                        unlist(str_split(cancer_gene, '\\|')),
                        unlist(str_split(ROI_hits, '\\|'))
                    )
                ) %>%
                mutate(check_score = check_scores$any_other_gene) %>%
                pull(check_score) %>% sum()
	  ) %>%
	  ungroup() 
}

annotate_precision.estimates <- function(tb, size_categories, precision_estimates) {
	tb %>%
		rowwise() %>%
		mutate(
			size_category = case_when(
				width >= size_categories$extreme.loh & CNV_type %!in% c('gain', 'loss') ~ 'extreme',
				width >= size_categories$extreme.cnv & CNV_type %in% c('gain', 'loss') ~ 'extreme',
				width >= size_categories$very.large.loh & CNV_type %!in% c('gain', 'loss') ~ 'very_large',
				width >= size_categories$very.large.cnv & CNV_type %in% c('gain', 'loss') ~ 'very_large',
				width >= size_categories$large.loh & CNV_type %!in% c('gain', 'loss') ~ 'large',
				width >= size_categories$large.cnv & CNV_type %in% c('gain', 'loss') ~ 'large',
				width >= size_categories$medium.loh & CNV_type %!in% c('gain', 'loss') ~ 'medium',
				width >= size_categories$medium.cnv & CNV_type %in% c('gain', 'loss') ~ 'medium',
				TRUE ~ 'small'
		  	),
			Precision_Estimate =
			  ifelse(CNV_type %in% c('gain', 'loss'),
				precision_estimates[[
					ifelse(n_initial_calls > 1, 'multiple_Callers', CNV_caller)]][[
					size_category]] +
				  (precision_estimates$Call_has_Gap * (probe_coverage_gap)) +
				  (precision_estimates$HighSNPDensity * (high_probe_density)),
				NA_real_
			  )
		) %>%
		ungroup() %>%
	  dplyr::select(-size_category)
}

annotate_call.label <- function(gr.or.tb, call_cat_config) {
    
    check_score.critical <- ifelse(is.null(call_cat_config$check_score.critical), NA, call_cat_config$check_score.critical)
    critical_excl <- call_cat_config$filters.exclude.critical 
    check_score.reportable <- ifelse(is.null(call_cat_config$check_score.reportable), NA, call_cat_config$check_score.reportable)
    check_score.reportable <- ifelse(is.na(check_score.reportable), check_score.critical, call_cat_config$check_score.reportable)
    reportable_excl <- call_cat_config$filters.exclude.reportable
    
    gr.or.tb %>%
        as_tibble() %>%
        mutate(
            Call_label = pmap_chr(
                list(reference_coverage, Check_Score, FILTER),
                \(ref_cov, check_score, FILTER) {
                    filters <- str_split(FILTER, ';') %>% unlist()
                    case_when(
                        !is.na(ref_cov)                        ~ 'Reference genotype',
                        check_score >= check_score.critical & 
                            !any(filters %in% critical_excl)   ~ 'Critical',
                        check_score >= check_score.reportable & 
                            !any(filters %in% reportable_excl) ~ 'Reportable',
                        TRUE                         		   ~ NA_character_
                    )
                }
            )
        )
}