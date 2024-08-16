suppressMessages(require(tidyverse))
suppressMessages(require(plyranges))

annotate_gaps <- function(gr, gapfile, min.perc.gap_area, gap_area.uniq_probes.rel, target_chrom_style = 'UCSC') {
	message('Annotation calls with gaps')

	gap_areas <- read_bed(gapfile) %>%
		select(-name, -score) %>%
        fix_CHROM_format(target_chrom_style) %>%
		mutate(gap_size = width) #%>%
		#filter(seqnames %in% get_chromosome_set())
	
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
	
	gr$percent_gap_coverage <- gap_ovs$perc_gap

	gap_slope <- gap_area.uniq_probes.rel[[1]]
	gap_intercept <- gap_area.uniq_probes.rel[[2]]

	gr <- gr %>%
		mutate(probe_coverage_gap =  ifelse(is.na(percent_gap_coverage),
										 FALSE,
										 percent_gap_coverage > min.perc.gap_area &
												(gap_slope * percent_gap_coverage + gap_intercept) <= log2(n_uniq_probes))
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

	gr$high_probe_density <- call_densities > density_value_cutoff

	gr
}
