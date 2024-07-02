combine_CNV_callers <- function(gr, min.greatest.region.overlap = 50, min.median.tool.coverage = 60) {
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
			reduce_ranges(ID = ID,
			              caller_merging_state = ifelse(plyranges::n() > 1, 'combined', 'no-overlap')) %>%
			as_tibble() %>%
			mutate(group.ID = paste('combined', CNV_type, seqnames, start, end, sep='_')) %>%
			dplyr::select(-strand, -seqnames) %>%
			rename_with(~paste0('group.', .), 1:3) %>%
			unnest(ID) %>%
			merge(as_tibble(gr), by = c('sample_id', 'CNV_type', 'ID')) %>%
			group_by(sample_id, CNV_type, pick(starts_with('group'))) %>%
			# The overlap here is *not* reciprocal overlap, as we need to deal with the possibility of one tool
			# 'merging' calls from another (or having >2 tools)
			mutate(overlap_merged_call = ifelse(caller_merging_state == 'combined', width / group.width * 100, NA),
			       max.ov = max(overlap_merged_call)) %>%
			group_by(sample_id, CNV_type, pick(starts_with('group')), CNV_caller) %>%
			# sum doesn't work if we have the filters in there and only group by CNV_caller
			mutate(tool.cov.sum = sum(overlap_merged_call)) %>%
			group_by(sample_id, CNV_type, seqnames, pick(starts_with('group'))) %>%
			mutate(tool.cov.ov.median = map2(list(CNV_caller), list(tool.cov.sum),
												\(CNV_caller, csum) tibble(CNV_caller = CNV_caller, csum = csum) %>% unique() %>%
																pull(csum) %>% median()
												),
			) %>% 
			# This makes testing *much* easier
			arrange(CNV_type, CNV_caller, start)

		#Now filter based on overlap of tool(CNV_caller) calls with merged region to check if the overlap can be accepted
		# - require at least 50% of the merged region to be covered by a single call to prevent chained overlaps
		# - require that the combined regions coverage from tools has a median of at least 60% (avg for only 2 tools)
		gr.changed <- ov %>%
			filter(caller_merging_state == 'combined' &
					   max.ov >= min.greatest.region.overlap &
					   tool.cov.ov.median >= min.median.tool.coverage
			) %>%
			select(-start, -end, -width,-ID) %>%
			rename_with(~str_remove(., '^group\\.')) %>%
			summarise(across(any_of(list_cols), ~list(.)),
					  overlap_merged_call = NA,
					  caller_merging_coverage = map2(CNV_caller, list(tool.cov.sum),
													\(CNV_caller, csum) tibble(t = CNV_caller, csum = csum) %>% unique() %>%
																	mutate(desc = paste0(t,'-', round(csum, 2))) %>%
																	pull(desc) %>% paste(collapse = ',')
													) %>% unlist(),
					  caller_merging_state = 'combined'
			) %>%
			as_granges() %>%
			ensure_list_cols()

		gr.unchanged <- ov %>%
			filter(caller_merging_state == 'no-overlap' |
					   max.ov < min.greatest.region.overlap |
					   tool.cov.ov.median < min.median.tool.coverage
			) %>%
			ungroup() %>%
			select(-starts_with('group.'), -max.ov, -tool.cov.sum, -tool.cov.ov.median) %>%
			mutate(caller_merging_coverage = NA_character_,
				   overlap_merged_call = NA,
			       caller_merging_state = 'no-overlap') %>%
			as_granges() %>%
			ensure_list_cols()

		gr.pre.overlap <-ov %>%
			filter(caller_merging_state == 'combined' &
					   max.ov >= min.greatest.region.overlap &
					   tool.cov.ov.median >= min.median.tool.coverage
			) %>%
			ungroup() %>%
			select(-starts_with('group.'), -max.ov, -tool.cov.sum, -tool.cov.ov.median) %>%
			mutate(caller_merging_coverage = NA_character_,
				   caller_merging_state = 'pre-overlap') %>%
			as_granges() %>%
			ensure_list_cols()

		return(bind_ranges(gr.changed, gr.unchanged, gr.pre.overlap))

	} else {
		gr <- gr %>%
			mutate(caller_merging_state = 'no-overlap',
				   #widths = NA,
				   overlap_merged_call = NA,
				   caller_merging_coverage = NA) %>%
			ensure_list_cols()

		return(gr)
	}

}