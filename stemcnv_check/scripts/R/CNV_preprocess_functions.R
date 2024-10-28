suppressMessages(require(tidyverse))
suppressMessages(require(plyranges))


merge_calls <- function(df.or.GR, merge_config, snp_vcf_gr) {
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
    
    rel.distance <- round(width(df.or.GR) * merge_config$call.extension.percent / 100)
    abs.distance <- rep(ceiling(merge_config$merge.gap.absolute/2), length(df.or.GR))
    
    if (merge_config$merge.gap.snps != 0) {
        # need to extend each call by i.e. 6 snps to close a gap of 10 snps (calls need to overlap for mrging) 
        probe_extend <- round(merge_config$merge.gap.snps/2)+1
        # get first & last SNP per CNV
        snps_by_cnv <- snp_vcf_gr %>%
            select(-ID) %>%
            mutate(i = 1:length(.)) %>%
            join_overlap_left(df.or.GR %>% select(ID)) %>%
            filter(!is.na(ID)) %>%
            group_by(ID)
        # Note: this will fail if any call has no probes (which should not occur with real data)
        if (length(unique(snps_by_cnv$ID)) != length(df.or.GR)) {
            stop('Some CNV calls have no SNP probes!')
        }        
        # extended from first/last snp per cnv by {probe_extend} SNPs
        # then get distance from that SNP to CNV start/end
        snp.dist.start <- snps_by_cnv %>% 
            filter(i==min(i)) %>% 
            mutate(
                ext_snp_i = i - probe_extend,
                ext_snp_i = ifelse(ext_snp_i < 1, 1, ext_snp_i),
                ext_snp_pos = snp_vcf_gr[ext_snp_i,] %>% start,
                snp.dist.start = start - ext_snp_pos,
                ID = factor(ID, levels = df.or.GR$ID)
            ) %>%
            # Need to extract distances in same order as input (df.or.GR)
            ungroup() %>%
            arrange(ID) %>%
            .$snp.dist.start
        snp.dist.end <- snps_by_cnv %>%
            filter(i==max(i)) %>% 
            mutate(
                ext_snp_i = i + probe_extend,
                ext_snp_i = ifelse(ext_snp_i > length(snp_vcf_gr), length(snp_vcf_gr), ext_snp_i),
                ext_snp_pos = snp_vcf_gr[ext_snp_i,] %>% end,
                snp.dist.end = ext_snp_pos - end ,
                ID = factor(ID, levels = df.or.GR$ID)
            ) %>%
            # Need to extract distances in same order as input (df.or.GR)
            ungroup() %>%
            arrange(ID) %>%
            .$snp.dist.end
        
        if(!all((start(df.or.GR) - snp.dist.start) %in% start(snp_vcf_gr))) {
            stop("SNP shifted start positions don't match existing SNP! Some CNV calls likely don't start at SNP probes.")
        }
        if(!all((end(df.or.GR) + snp.dist.end) %in% start(snp_vcf_gr))) {
            stop("SNP shifted end positions don't match existing SNP! Some CNV calls likely don't end at SNP probes.")
        }    
    } else {
        snp.dist.start <- 0
        snp.dist.end <- 0
    }
    # snp distance will differ for start & end, so need to calculate separately
    df.or.GR$start_reduce <- map2_int(
        pmap_int(list(rel.distance, abs.distance, snp.dist.start), max),
        merge_config$maximum.gap.allowed, 
        min, na.rm = TRUE
    )
    df.or.GR$end_increase <- map2_int(
        pmap_int(list(rel.distance, abs.distance, snp.dist.end), max),
        merge_config$maximum.gap.allowed, 
        min, na.rm = TRUE
    )

    df.or.GR %>%
        mutate(
            orig_start = start,
            orig_end = end,
            start = start - start_reduce,
            end = end + end_increase,
        ) %>%
        group_by(sample_id, CNV_caller, CNV_type) %>%
        mutate(short_ID = paste(str_extract(ID, '[0-9]+_[0-9]+$'), CN, sep='_CN')) %>%
        reduce_ranges(
            n_initial_calls = plyranges::n(),
            # Collect initial start, end & CN
            initial_call_details = ifelse(
                plyranges::n() > 1,
                paste(short_ID, collapse = ';'),
                NA_character_
            ),
            CN = median(CN),
            orig_start = orig_start,
            orig_end = orig_end
        ) %>%
        # Restore the original start & end coords, but ONLY from the ends of new merged calls
        mutate(
            start = min(orig_start), 
            end = max(orig_end)
        ) %>%
        select(-orig_start, -orig_end) %>%
        mutate(ID = paste(CNV_caller, CNV_type, seqnames, start, end, sep='_')) %>%
        add_snp_probe_counts(snp_vcf_gr)
}

# n_nsp & uniq_snp annotation
add_snp_probe_counts <- function(gr, snp_vcf_gr) {
    message('adding number of SNP probes to calls')
    
    # Return early if input is empty
    if (length(gr) == 0) {
        return(
            gr %>%
                mutate(
                    n_probes = integer(),
                    n_uniq_probes = integer(),
                    probe_density_Mb = double(),
                )
        )
    }
    
    if (length(snp_vcf_gr) == 0) {
        # This is for testing
        warning('No SNP probes provided for annotation!')
    } else if (seqlevelsStyle(gr) %>% head(1) != seqlevelsStyle(snp_vcf_gr) %>% head(1)) {
        stop('seqlevelsStyle of input (CNV) GRanges and SNP ref do not match!')
    }
    
    # Note: with default filter settings, these two WILL be the same
	gr$n_probes <- count_overlaps(gr, snp_vcf_gr)
	gr$n_uniq_probes <- count_overlaps(
        gr,
        # reduce_ranges without grouping will also merge probes NEXT to each other 
        snp_vcf_gr %>%
            group_by(start) %>% 
            reduce_ranges()
    )
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
            f_size = ifelse(width < min.length, 'min_size', NA_character_),
            f_probes = ifelse(n_probes < min.snp, 'min_probes', NA_character_),
            f_density = ifelse(probe_density_Mb < min.snp.density, 'min_density', NA_character_),
            FILTER = ifelse(
                is.na(f_size) & is.na(f_probes) & is.na(f_density),
                NA_character_,
                paste(na.omit(c(f_size, f_probes, f_density)), collapse = ';'))
        ) %>%
        select(-f_size, -f_probes, -f_density) %>%
        ungroup() %>%
        as_granges()
}

apply_preprocessing <- function(cnv_gr, snp_vcf_gr, tool_config) {
        
    cnv_gr %>%
        # Properly sort by updated seqnames
        sort() %>%
        # Combined nearby calls & update SNP counts
        merge_calls(tool_config$call.merging, snp_vcf_gr) %>%
        # Add FILTER column
        add_call_prefilters(tool_config) %>%
        sort()
    
}

# FIXME (future): some way to get (PennCNV) GC correction for LRR?
# FIXME (future): add solution for BAF
get_median_LRR <- function(gr, snp_vcf_gr) {
    
    median_lrr <- snp_vcf_gr %>%
        plyranges::select(LRR) %>%
        join_overlap_inner(gr) %>%        
        as_tibble() %>%
        group_by(ID) %>%
        summarise(LRR = median(LRR) %>% round(3))
    
    # If calls overlap (i.e. not meeting merging criteria) using plyranges here would be problematic
    as_tibble(gr) %>%
        full_join(median_lrr, by = 'ID') %>%
        as_granges()
}