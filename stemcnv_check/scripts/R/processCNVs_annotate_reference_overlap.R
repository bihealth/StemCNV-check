suppressMessages(require(tidyverse))
suppressMessages(require(plyranges))

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
