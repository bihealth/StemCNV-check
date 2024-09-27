suppressMessages(require(tidyverse))
suppressMessages(library(GenomicRanges))
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
