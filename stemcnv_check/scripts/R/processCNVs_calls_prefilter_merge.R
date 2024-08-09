merge_calls <- function(df.or.GR, merge.distance, snp.vcf) {
    message('Merging nearby raw calls')
    if (is.data.frame(df.or.GR)) {
        df.or.GR <- df.or.GR %>%
		    as_granges()
	}
    # Return early if input is empty
    if (length(df.or.GR) == 0) {
        return(
            df.or.GR %>%
                # ensure same order as with actual data
                select(-CN, -ID) %>%
                mutate(
                    n_initial_calls = integer(),
                    initial_call_details = character(),
                    CN = integer(),
                    ID = character(),
                    n_probes = integer(),
                    n_uniq_probes = integer(),
                    probe_density_Mb = double(),
                )
        )
    }
    
	df.or.GR %>%
		group_by(sample_id, CNV_caller, CNV_type) %>%
		stretch(merge.distance) %>%
        mutate(short_ID = paste(str_extract(ID, '[0-9]+_[0-9]+$'), CN, sep='_CN')) %>%
		reduce_ranges(
			n_initial_calls = plyranges::n(),
            # Collect initial start, end & CN
			initial_call_details = ifelse(
                plyranges::n() > 1,
                paste(short_ID, collapse = ','),
                NA_character_
            ),
			CN = median(CN),
			# Too unreliable to actuallly use
            # caller_confidence = median(caller_confidence)
			) %>%
		stretch(-1*merge.distance) %>%
		mutate(ID = paste(CNV_caller, CNV_type, seqnames, start, end, sep='_')) %>%
        add_snp_probe_counts(snp.vcf)
}

# n_nsp & uniq_snp annotation
add_snp_probe_counts <- function(gr, snp.vcf) {
    message('adding number of SNP probes to calls')

    # VCF POS should be 1-based,
    # Granges are also 1-based
    # AND are (fully) open [= start & end are included]
    snp_probes_gr <- vcfR2tidy(snp.vcf, info_only = T) %>%
        .$fix %>%
        as_granges(seqnames = CHROM, start = POS, width = 1) %>%
        plyranges::filter(FILTER == 'PASS') %>%
        # reduce_ranges without grouping will also merge probes NEXT to each other 
        group_by(start)
    
    # Note: with default filter settings, these two WILL be the same
	gr$n_probes <- count_overlaps(gr, snp_probes_gr)
    
	gr$n_uniq_probes <- count_overlaps(gr, reduce_ranges(snp_probes_gr))
    gr$probe_density_Mb <- gr$n_uniq_probes / width(gr) * 1e6

	gr

}


add_call_prefilters <- function(gr, tool_config) {
	message('Pre-filtering calls')
    
    # Return early if input is empty
    if (length(gr) == 0) {
        return(
            gr %>%
                mutate(
                    FILTER = character(),
                )
        )
    }
    
    min.snp <- tool_config$filter.minprobes
    min.length <- tool_config$filter.minlength
    min.snp.density <- tool_config$filter.mindensity.Mb
   
    gr %>% 
        as_tibble() %>%
        rowwise() %>%
        mutate(
            f_size = ifelse(width < min.length, 'Size', NA_character_),
            f_probes = ifelse(n_probes < min.snp, 'n_probes', NA_character_),
            f_density = ifelse(probe_density_Mb < min.snp.density, 'Density', NA_character_),
            FILTER = ifelse(
                is.na(f_size) & is.na(f_probes) & is.na(f_density),
                'PASS',
                paste(na.omit(c(f_size, f_probes, f_density)), collapse = ';'))
        ) %>%
        select(-f_size, -f_probes, -f_density) %>%
        ungroup() %>%
        as_granges()
}


