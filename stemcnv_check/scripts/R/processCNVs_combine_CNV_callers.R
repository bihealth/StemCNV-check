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
			    caller_merging_state = ifelse(plyranges::n() > 1, 'combined', 'no-overlap')
            ) %>%
			as_tibble() %>%
			mutate(group.ID = paste('StemCNV-check', CNV_type, seqnames, start, end, sep='_')) %>%
			dplyr::select(-strand, -seqnames) %>%
			rename_with(~paste0('group.', .), 1:3) %>%
			unnest(ID) %>%
			merge(as_tibble(gr), by = c('sample_id', 'CNV_type', 'ID')) %>%
			group_by(sample_id, CNV_type, pick(starts_with('group'))) %>%
			# The overlap here is *not* reciprocal overlap, as we need to deal with the possibility of one tool
			# 'merging' calls from another (or having >2 tools)
			mutate(
                overlap_merged_call = ifelse(caller_merging_state == 'combined', width / group.width * 100, NA_real_),
			    max.ov = max(overlap_merged_call)
            ) %>%
			group_by(sample_id, CNV_type, pick(starts_with('group')), CNV_caller) %>%
			# sum doesn't work if we have the filters in there and only group by CNV_caller
			mutate(tool.cov.sum = sum(overlap_merged_call)) %>%
			group_by(sample_id, CNV_type, seqnames, pick(starts_with('group'))) %>%
			mutate(
                tool.cov.ov.median = map2(
                    list(CNV_caller), list(tool.cov.sum),
					\(CNV_caller, csum) 
                        tibble(CNV_caller = CNV_caller, csum = csum) %>% 
                            unique() %>% pull(csum) %>% median()
                ),
			) %>% 
			# This makes testing *much* easier
			arrange(CNV_type, CNV_caller, start)

		#Now filter based on overlap of tool(CNV_caller) calls with merged region to check if the overlap can be accepted
		# - require at least 50% of the merged region to be covered by a single call to prevent chained overlaps
		# - require that the combined regions coverage from tools has a median of at least 60% (avg for only 2 tools)
		gr.changed <- ov %>%
			filter(
                caller_merging_state == 'combined' &
				max.ov >= processing_config$tool.overlap.greatest.call.min.perc &
				tool.cov.ov.median >= processing_config$tool.overlap.median.cov.perc
			) %>%
			# select(-start, -end, -width, -ID) %>%
			summarise(
                # Add new columns
				caller_merging_coverage = map2(
                    list(CNV_caller),
                    list(tool.cov.sum),
                    \(CNV_caller, csum) 
                        tibble(t = CNV_caller, csum = csum) %>% unique() %>%
                            mutate(desc = paste0(t,'-', round(csum, 2))) %>%
                            pull(desc) %>% paste(collapse = ';')
                ) %>% unlist(),
				caller_merging_state = 'combined',
                #FIXME: in the long run this can probably replace keeping the changed calls 
                # and remove the "caller_merging_state" column
                n_initial_calls = dplyr::n(),
                # Collect initial start, end, CN, callers & coverages
                #FIXME: maybe this should be a list column?
                initial_call_details = paste(
                    str_glue("{CNV_caller}_{start}-{end}_CN{CN}_cov{round(overlap_merged_call, 2)}"), 
                    collapse = '|'
                ),
                # Transform existing columns
                across(any_of(get_list_cols()), ~list(.)),
                CN = median(CN),
                overlap_merged_call = NA_real_,
			) %>%
            rename_with(~str_remove(., '^group\\.')) %>%
			as_granges() %>%
            # this should not be necessary?
			ensure_list_cols() %>%
            add_snp_probe_counts(snp_vcf)

		gr.unchanged <- ov %>%
			filter(
                caller_merging_state == 'no-overlap' |
				max.ov < processing_config$tool.overlap.greatest.call.min.perc |
				tool.cov.ov.median < processing_config$tool.overlap.median.cov.perc
			) %>%
			ungroup() %>%
			select(-starts_with('group.'), -max.ov, -tool.cov.sum, -tool.cov.ov.median) %>%
			mutate(
                overlap_merged_call = NA_real_,
				caller_merging_coverage = NA_character_,
			    caller_merging_state = 'no-overlap',
                # From here on this is relative to stemcnv-check tool, not orig. tool
                initial_call_details = NA_character_,
                n_initial_calls = 1
            ) %>%
			as_granges() %>%
			ensure_list_cols()

		gr.pre.overlap <-ov %>%
			filter(
                caller_merging_state == 'combined' &
				max.ov >= processing_config$tool.overlap.greatest.call.min.perc &
				tool.cov.ov.median >= processing_config$tool.overlap.median.cov.perc
			) %>%
			ungroup() %>%
			select(-starts_with('group.'), -max.ov, -tool.cov.sum, -tool.cov.ov.median) %>%
			mutate(
                caller_merging_coverage = NA_character_,
				caller_merging_state = 'pre-overlap',
                # From here on this is relative to stemcnv-check tool, not orig. tool
                initial_call_details = NA_character_,
                n_initial_calls = NA_integer_,
            ) %>%
			as_granges() %>%
			ensure_list_cols()

		return(bind_ranges(gr.changed, gr.unchanged, gr.pre.overlap))

	} else {
		gr <- gr %>%
			#plyranges::mutate fails on empty gr
			as_tibble() %>%
			mutate(
				overlap_merged_call = NA_real_,
				caller_merging_coverage = NA_character_,
				caller_merging_state = 'no-overlap',
                # From here on this is relative to stemcnv-check tool, not orig. tool
                initial_call_details = NA_character_,
                n_initial_calls = 1
                
			) %>%
			ensure_list_cols()

		return(gr)
	}

}