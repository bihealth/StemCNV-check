# n_nsp & uniq_snp annotation
get_accurate_snp_probe_count <- function(gr, unfiltered.snps.file) {
	message('(re)calculating number of SNP probes per call')

	snp_probes_gr <- read_raw(unfiltered.snps.file) %>%
		as_granges(seqnames = Chr, start = Position, width = 1)

	n_snp <- count_overlaps(gr, snp_probes_gr)
	n_snp_uniq <- count_overlaps(gr, reduce_ranges(snp_probes_gr))

	gr$n_snp_probes <- n_snp
	gr$n_uniq_probe_positions <- n_snp_uniq

	gr

}

annotate_gaps <- function(gr, gapfile, min.perc.gap_area, gap_area.uniq_probes.rel) {
	message('Annotation calls with gaps')

	gap_areas <- read_bed(gapfile) %>%
		select(-name, -score) %>%
		mutate(gap_size = width) %>%
		filter(seqnames %in% get_chromosome_set())
	
	gap_ovs <- join_overlap_left(gr, gap_areas) %>%
		group_by(ID) %>%
		reduce_ranges(gap_size_sum = sum(gap_size)) %>%
		mutate(perc_gap = ifelse(is.na(gap_size_sum), 0, gap_size_sum  / width)) 
	
	# Ensure that calls _fully_ overlap gaps
	# This should be given since gaps are defined as / derived from large inter-probe distances
	if (find_overlaps_within(gap_areas, gap_ovs) %>% length() != length(gap_ovs)) {
		stop(str_glue('Gaps defined by "{gapfile}" are not contained within CNV calls!'))
	}
	if (find_overlaps_within(gr, gap_areas) %>% length() > 0) {
		stop(str_glue('Gaps defined by "{gapfile}" fully overlap some CNV calls!'))
	}
	
	gr$percent_gap_coverage <- gap_ovs$perc_gap

	gap_slope <- gap_area.uniq_probes.rel[[1]]
	gap_intercept <- gap_area.uniq_probes.rel[[2]]

	gr <- gr %>%
		mutate(probe_coverage_gap =  ifelse(is.na(percent_gap_coverage),
										 FALSE,
										 percent_gap_coverage > min.perc.gap_area &
												(gap_slope * percent_gap_coverage + gap_intercept) <= log2(n_uniq_probe_positions))
		)

	return(gr)
}

annotate_high_density <- function(gr, density_file, density.quantile.cutoff) {
	message('Annotation calls for high probe density')

	array_density <- read_bed(density_file) %>%
		select(-name) %>%
		mutate(density = score) %>%
		filter(seqnames %in% get_chromosome_set()) %>%
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
