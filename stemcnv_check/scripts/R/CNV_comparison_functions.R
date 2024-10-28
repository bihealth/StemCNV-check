suppressMessages(require(tidyverse))
suppressMessages(require(plyranges))

load_cnv_callers <- function(input_vcf_files) {
    lapply(input_vcf_files, parse_cnv_vcf) %>%
        bind_ranges()
}


combine_CNV_callers <- function(gr, processing_config, snp_vcf) {
	message('Combining calls from multiple callers')

	ov_test <- gr %>%
		group_by(sample_id, CNV_type) %>%
		reduce_ranges()

	# If overlaps exist (if not reduce_ranges can't make proper list cols):
	# Merge any number of overlaps for the same CNV.state together, summarise metadata into list cols
	# & calculate coverage of the final merged call by individual ones
	if (length(ov_test) < length(gr))  {
		ov <- gr %>%
			group_by(sample_id, CNV_type) %>%
			reduce_ranges(
                ID = ID,
                n_initial_calls = plyranges::n(),
            ) %>%
			as_tibble() %>%
			mutate(group.ID = paste('StemCNV-check', CNV_type, seqnames, start, end, sep='_')) %>%
			dplyr::select(-strand) %>%
			rename_with(~paste0('group.', .), 1:4) %>%
			unnest(ID) %>%
			full_join(
                as_tibble(gr) %>% select(-n_initial_calls), 
                by = c('sample_id', 'CNV_type', 'ID')
            ) %>%
			group_by(sample_id, CNV_type, pick(starts_with('group'))) %>%
			# Get the overlap of each call to the complete area, then get 
            # - the max ov per area from a single call
            # - the minimum of the summed area coverage for each caller
            # (recripocal coverage does not apply, because >2 calls (or callers) might be involved)
			mutate(
                overlap_merged_call = width / group.width * 100,
			    max.ov = max(overlap_merged_call),
                min.tool.cov.sum = map2(
                    list(CNV_caller),
                    list(overlap_merged_call),
                    \(CNV_caller, ov) {
                        tibble(caller = CNV_caller, ov = ov) %>%
                            group_by(caller) %>%
                            summarise(sum = sum(ov)) %>%
                            pull(sum) %>% min()
                    }
                )
            ) %>%
			# This makes testing *much* easier
			arrange(CNV_type, CNV_caller, start)

		#Now filter based on overlap of tool(CNV_caller) calls with merged region to check if the overlap can be accepted
		# - require at least 50% of the merged region to be covered by a single call to prevent chained overlaps
		# - require that the combined regions coverage from each tool has a sum of at least 60%
		gr.changed <- ov %>%
			filter(
                n_initial_calls > 1 &
				max.ov >= processing_config$tool.overlap.greatest.call.min.perc &
				min.tool.cov.sum >= processing_config$tool.overlap.min.cov.sum.perc
			) %>%
            # Add new columns
			summarise(
                n_initial_calls = dplyr::n(),
                # Collect initial callers, start, end, CN, FILTER & coverages
                initial_call_details = paste(
                    str_glue(
                        "{CNV_caller}_{start}-{end}_CN{CN}_",
                        "cov{round(overlap_merged_call, 2)}_",
                        "{ifelse(is.na(FILTER), 'PASS', str_replace_all(FILTER, ';', '&'))}"
                    ), 
                    collapse = '|'
                ),
                CNV_caller = 'StemCNV-check',
                CN = median(CN),
			) %>%
            rename_with(~str_remove(., '^group\\.')) %>%
			as_granges() %>%
            add_snp_probe_counts(snp_vcf) %>%
            add_call_prefilters(processing_config)

		gr.unchanged <- ov %>%
			filter(
                n_initial_calls == 1 |
				max.ov < processing_config$tool.overlap.greatest.call.min.perc |
				min.tool.cov.sum < processing_config$tool.overlap.min.cov.sum.perc
			) %>%
			ungroup() %>%
			select(
                -starts_with('group.'),
                -max.ov, -min.tool.cov.sum, -overlap_merged_call
            ) %>%
			mutate(
                # From here on this is relative to stemcnv-check tool, not orig. tool
                initial_call_details = NA_character_,
                n_initial_calls = 1
            ) %>%
			as_granges() 

		return(bind_ranges(gr.changed, gr.unchanged))

	} else {
		gr <- gr %>%
			#plyranges::mutate fails on empty gr
			as_tibble() %>%
			mutate(
                # From here on this is relative to stemcnv-check tool, not orig. tool
                initial_call_details = NA_character_,
                n_initial_calls = 1
			) %>%
            as_granges()

		return(gr)
	}

}


annotate_reference_overlap <- function(gr_in, gr_ref, min.reciprocal.coverage.with.ref = 0.8) {
	message('Annotation calls with reference overlap')
	#If pair_overlaps is an empty set/df the min/max functions will give a warning, while later Granges functions would crash
	# dplyr::summarise will also give a (deprecation) warning if it gets a fully empty tibble
	suppressWarnings(
	gr <- pair_overlaps(gr_in, gr_ref) %>%
			as_tibble() %>% rowwise() %>%
			mutate(
                width.combined = max(granges.x.end, granges.y.end) - min(granges.x.start, granges.y.start) + 1,
				width.overlap = min(granges.x.end, granges.y.end) - max(granges.x.start, granges.y.start) + 1,
				reference_coverage = round(100 * width.overlap / granges.x.width, 2),
				cov.ref.by.sample = round(100 * width.overlap / granges.y.width, 2),
				# need to flatten/fully unpack tool.y or we get list of lists
				reference_caller = list(unlist(CNV_caller.y)),
				ref.state = CNV_type.y
			) %>%
			dplyr::select(-contains('.y'), -contains('width')) %>%
			filter(reference_coverage >= min.reciprocal.coverage.with.ref &
					   cov.ref.by.sample >= min.reciprocal.coverage.with.ref &
					   CNV_type.x == ref.state)
	)

	if (nrow(gr) > 0) {
		gr <- gr %>%
			group_by(granges.x.seqnames, granges.x.start, granges.x.end, CNV_caller.x, CNV_type.x) %>%
			summarise(
                across(ends_with('.x'), ~ unique(.)),
				reference_overlap = TRUE,
				reference_coverage = sum(reference_coverage),
				reference_caller = paste(unique(unlist(reference_caller)), collapse = ';'),
			) %>%
			dplyr::rename_with(~ str_remove(., '.x$') %>% str_remove('^granges.x.')) %>%
			# Preserve original column order
			select(any_of(colnames(as_tibble(gr_in))), 
						 reference_overlap, reference_coverage, reference_caller
			) %>%
			as_granges()

		# `gr` only has (filtered) overlaps, need to rebuild the full callset if matching ref calls were found
		# Get Original regions not in the merged call set by ID
		non.ovs <- gr_in %>% filter(ID %!in% gr$ID)
		# mutate fails on empty GRanges object
		if (length(non.ovs) > 0) { non.ovs <- non.ovs %>% 
			mutate(
                reference_overlap = FALSE,
			    reference_coverage = NA_real_,
				reference_caller = NA_character_
			)
		}

		gr_out <- bind_ranges(
			# Calls with matching reference
			gr,
			# Remaining calls not (sufficiently) overlapping reference
			non.ovs
		)

	} else {
		gr_out <- gr_in %>%
			as_tibble() %>%
			mutate(
                reference_overlap = FALSE,
				reference_coverage = NA_real_,
				reference_caller = NA_character_
            ) %>%
			as_granges()
	}
	
	gr_out

}
