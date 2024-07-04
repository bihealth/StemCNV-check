merge_calls <- function(df.or.GR, merge.distance) {
	message('Merging nearby raw calls')
	if (is.data.frame(df.or.GR)) {
		df.or.GR <- df.or.GR %>%
			dplyr::rename(seqnames = Chr) %>%
			as_granges()
	}
	df.or.GR %>%
		group_by(CNV_type, sample_id, CNV_caller) %>%
		stretch(merge.distance) %>%
		reduce_ranges(
			n_premerged_calls = plyranges::n(),
			n_snp_probes = sum(n_snp_probes),
			copynumber = paste(unique(copynumber), collapse = ','),
			caller_confidence = median(caller_confidence)
			) %>%
		stretch(-1*merge.distance) %>%
		mutate(ID = paste(CNV_caller, CNV_type, seqnames, start, end, sep='_'))
}

prefilter_calls <- function(df.or.GR, min.snp, min.length, min.snp.density) {
	message('Pre-filtering calls')
	if (is.data.frame(df.or.GR)){
		df.or.GR <- dplyr::filter(df.or.GR, n_snp_probes >= min.snp & length >= min.length & snp.density >= min.snp.density)
	} else {
		df.or.GR <- plyranges::filter(df.or.GR,	n_snp_probes >= min.snp & width >= min.length & (n_snp_probes / width * 1e6) >= min.snp.density )
	}
	df.or.GR
}


