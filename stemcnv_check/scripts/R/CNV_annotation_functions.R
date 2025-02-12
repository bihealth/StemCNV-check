suppressMessages(require(tidyverse))
suppressMessages(library(GenomicRanges))
suppressMessages(require(plyranges))

annotate_impact_lists <- function(gr, hotspot_gr, list_name) {
	message('Annotation calls with gene lists')

	# Make an extra col listing all overlapping directly defined genes
	gr@elementMetadata[[list_name]] <- NA_character_
	ov <- group_by_overlaps(gr, hotspot_gr)
	if (length(ov) > 0) {
		ov_hits <- ov %>%
			mutate(type_check = str_detect(CNV_type, str_replace(call_type, '^any$', '.*'))) %>%
			filter(type_check) %>%
			reduce_ranges(hits = paste(unique(sort(hotspot)),collapse = '|'))
		gr[ov_hits$query,]@elementMetadata[[list_name]] <- ov_hits$hits
	}
	return(gr)
}

#FIXME: this could probably be replaced/refactored into annotate_impact_lists
annotate_roi <- function(gr, roi_tb, gr_genes, gr_info, config) {
	message('Annotating calls with ROI')

    target_chrom_style <- get_target_chrom_style(config, gr)
	roi_gr <- get_roi_gr(roi_tb, gr_genes, gr_info, target_chrom_style)

	gr$ROI_hits <- NA_character_
	ov <- group_by_overlaps(gr, roi_gr)
	if (length(ov) > 0) {
		ov_hits <- ov %>%
			reduce_ranges(hits = paste(unique(hotspot), collapse = '|'))
		gr[ov_hits$query,]$ROI_hits <- ov_hits$hits
	}
	gr
}


annotate_gaps <- function(gr, gapfile, min.perc.gap_area, gap_area.uniq_probes.rel, target_chrom_style = 'UCSC') {
	message('Annotation calls with gaps')
    
	gap_areas <- read_bed(gapfile) %>%
		select(-name, -score) %>%
        fix_CHROM_format(target_chrom_style) %>%
		mutate(gap_size = width)
    
	# Sanity check: Are any CNV breakpoints in gap?
	# This should never happen since CNVs & gaps are derived from probes
	start_or_end_in_gap <- bind_ranges(
		gr %>% anchor_start() %>% mutate(width = 1),
		gr %>% anchor_end() %>% mutate(width = 1)
	) %>% filter_by_overlaps(gap_areas)
	if ( length(start_or_end_in_gap) > 0 ) {
		call_ids <- start_or_end_in_gap$ID %>% unique() %>% paste(collapse = ', ')
		stop(str_glue('CNV call endpoint(s) overlap with gap areas from "{gapfile}": {call_ids}'))
	}
	
	gap_ovs <- join_overlap_left(gr, gap_areas) %>%
		group_by(ID) %>%
		reduce_ranges(gap_size_sum = sum(gap_size)) %>%
		mutate(perc_gap = ifelse(is.na(gap_size_sum), 0, gap_size_sum  / width))
	
	gr$Gap_percent <- gap_ovs$perc_gap

	gap_slope <- gap_area.uniq_probes.rel[[1]]
	gap_intercept <- gap_area.uniq_probes.rel[[2]]

	gr <- gr %>%
		mutate(
            probe_coverage_gap = ifelse(
                is.na(Gap_percent),
                FALSE,
                Gap_percent > min.perc.gap_area & 
                    (gap_slope * Gap_percent + gap_intercept) <= log2(n_uniq_probes)
            ),
            FILTER = map2_chr(FILTER, probe_coverage_gap, \(f, gaps) {
                gaps[is.na(gaps)] <- FALSE
                ifelse(
                    gaps,
                    paste(na.omit(c(f, 'probe_gap')), collapse = ';'),
                    f
                )
            })
		)
	return(gr)
}


annotate_high_density <- function(gr, density_file, density.quantile.cutoff, target_chrom_style = 'UCSC') {
	message('Annotation calls for high probe density')

	array_density <- read_bed(density_file) %>%
		select(-name) %>%
        fix_CHROM_format(target_chrom_style) %>%
		mutate(density = score) %>%
		# filter(seqnames %in% get_chromosome_set()) %>%
		filter(density > 0)

	density_value_cutoff <- quantile(array_density$density, density.quantile.cutoff)

	call_densities <- join_overlap_left(gr, array_density) %>%
		group_by(ID) %>%
		# partially overlapping density windows might over-proportionally affect calls (since we don't the density for % overlapped)
		# However, this only matters for very small calls (and density in neighbouring windows will likely be similar/correlated anyway)
		reduce_ranges(density = mean(density)) %>%
		as_tibble() %>%
		pull(density)

	gr %>%
        mutate(
            high_probe_density = call_densities > density_value_cutoff,
            FILTER = map2_chr(FILTER, high_probe_density, \(f, dens) {
                dens[is.na(dens)] <- FALSE
                ifelse(
                    dens,
                    paste(na.omit(c(f, 'high_probe_dens')), collapse = ';'),
                    f
                )
            })
        )
}


# Sanitize output & add gene overlap annotation
annotate_gene_overlaps <- function(gr, gr_genes) {

	gr$n_genes <- count_overlaps(gr, gr_genes)
	gr$overlapping_genes <- NA_character_
	ov_genes <- group_by_overlaps(gr, gr_genes) %>% reduce_ranges(genes = paste(gene_name, collapse = '|'))
	gr[ov_genes$query,]$overlapping_genes <- ov_genes$genes

	return(gr)
}


annotate_cnv.check.score <- function(tb, stemcell_hotspots_gr, dosage_sensitive_gene_gr, cancer_genes_gr, check_scores, sample_sex) {

    cn1_3 <- check_scores$single_copy_factor
    cn0_4 <- check_scores$double_copy_factor
    cn2   <- check_scores$neutral_copy_factor
    flat  <- check_scores$flat_decrease  
    
	tb %>% 
        rowwise() %>% 
        mutate(
            Check_Score =
                case_when(
                    # Male sex chromosomes can not have LOH, and should always be treated a single copy CNVs (unless double gain)
                    sample_sex == 'm' & str_detect(seqnames, 'X|Y')
                        & CN <= 2                                         ~ cn1_3 * log(width) * log(width) - flat,
                    sample_sex == 'm' & str_detect(seqnames, 'X|Y') 
                        & CN > 2                                          ~ cn0_4 * log(width) * log(width) - flat,
                    # large CN (0/4)
                    CNV_type %in% c('gain', 'loss') & (CN == 0 | CN == 4) ~ cn0_4 * log(width) * log(width) - flat,
                    # single gain/loss CN (1/3)
                    CNV_type %in% c('gain', 'loss')                       ~ cn1_3 * log(width) * log(width) - flat,
                    # LOH == CN2
                    CNV_type == 'LOH'                                     ~   cn2 * log(width) * log(width) - flat
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
    
    # local function that can refer to call_cat_config
    assign_single_call_label <- function(reference_overlap, check_score, FILTER) {
        # Split single string into list of filters
        filters <- str_split(FILTER, ';') %>% unlist()

        for (cat in names(call_cat_config)) {
            # Allow minimum_check_score to be unset (NULL), in which case it is set to -inf
            min_check_score <- ifelse(
                is.null(call_cat_config[[cat]]$minimum_check_score),
                -Inf,
                call_cat_config[[cat]]$minimum_check_score
            )
            # check if the call matches the category
            if (
                check_score >= min_check_score & 
                !any(filters %in% call_cat_config[[cat]]$not_allowed_vcf_filters) & 
                reference_overlap == call_cat_config[[cat]]$reference_match
            ) {
                return(cat)
            }
        }
        message(
            'No category found for call with reference_coverage: ', reference_overlap, ', Check_Score = ', check_score, 
            ', and FILTERs: ', FILTER, '. Falling back to last defined category: ', cat
        )
        return(cat)
    }
    
    gr.or.tb %>%
        as_tibble() %>%
        mutate(
            Call_label = pmap_chr(list(reference_overlap, Check_Score, FILTER), assign_single_call_label)
        )
}