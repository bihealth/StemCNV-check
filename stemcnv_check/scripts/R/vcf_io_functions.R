suppressMessages(require(tidyverse))
suppressMessages(require(plyranges))
suppressMessages(require(vcfR))

# VCF reading
vcfR_to_tibble <- function(vcf, info_fields = NULL, format_fields = NULL){
    #Note: this will fail if format fields are not properly annotated in header
    listobj <- vcfR2tidy(vcf, alleles = FALSE, gt_column_prepend = '', single_frame = TRUE,
                         info_fields = info_fields, format_fields = format_fields)
    listobj$dat %>% dplyr::rename(sample_id = Indiv)       
}

parse_snp_vcf <- function(vcf, 
                         info_fields = FALSE, 
                         format_fields = c("LRR", "BAF"), 
                         apply_filter = TRUE) {
    if (typeof(vcf) == 'character') {
        vcf <- read.vcfR(vcf, verbose = F)
    }
    # VCF POS should be 1-based,
    # Granges are also 1-based
    # AND are (fully) open [= start & end are included]
    vcf <- vcf %>%
        vcfR_to_tibble(info_fields = info_fields, format_fields = format_fields) %>%
        mutate(sample_id = str_remove(sample_id, '\\.gencall')) %>%
        as_granges(seqnames = CHROM, start = POS, width = 1)
    if (apply_filter) {
        vcf <- vcf %>% filter(FILTER == 'PASS')
    }
    return(vcf)
}

# parse_cnv_vcf
parse_cnv_vcf <- function(vcf, 
                         info_fields = NULL, 
                         format_fields = NULL, 
                         apply_filter = TRUE) {
    if (typeof(vcf) == 'character') {
        vcf <- read.vcfR(vcf, verbose = F)
    }
    # VCF POS should be 1-based,
    # Granges are also 1-based
    # AND are (fully) open [= start & end are included]
    vcf <- vcf %>%
        vcfR_to_tibble(info_fields = info_fields, format_fields = format_fields) %>%
        mutate(
            CNV_type = str_remove(ALT, '<') %>% str_remove('>') %>%
                str_replace('DUP', 'gain') %>%
                str_replace('DEL', 'loss') %>%
                str_replace('CNV:LOH', 'LOH'),
            CNV_caller = str_extract(TOOL, '(?<=caller=)[^;]+'),
            n_initial_calls = str_extract(TOOL, '(?<=n_initial_calls=)[^;]+'),
            initial_call_details = str_extract(TOOL, '(?<=initial_call_details=)[^;]+'),
        ) %>%
        rename_with(~str_to_lower(.), contains('PROBE')) %>%
        dplyr::rename(probe_density_Mb = probe_dens) %>%
        select(-REF, -ALT, -QUAL, -SVCLAIM, -TOOL) %>%
        as_granges(seqnames = CHROM, start = POS + 1, end = END, width = SVLEN)
    if (apply_filter) {
        vcf <- vcf %>% filter(FILTER == 'PASS')
    }
    return(vcf)
}


# CNV vcf writing
static_cnv_vcf_header <- function(toolconfig, extra_annotation = FALSE, INFO = TRUE, FORMAT = TRUE, FILTER = TRUE) {

    info <- c(
        '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the longest variant described in this record">',
	    '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of structural variant">',
        '##INFO=<ID=SVCLAIM,Number=A,Type=String,Description="Claim made by the structural variant call. Valid values are D, J, DJ for abundance, adjacency and both respectively">',
        # Maybe use Number=. for these, if we want to allow multiple values after call overlaps
        '##INFO=<ID=N_PROBES,Number=1,Type=Integer,Description="Number of array probes in segment">',
        '##INFO=<ID=N_UNIQ_PROBES,Number=1,Type=Integer,Description="Number of unique array probe positions in segment">',
        '##INFO=<ID=PROBE_DENS,Number=1,Type=Float,Description="Density of Probes in segment (Probes / 10Mb)">'
    )
    if (extra_annotation) {
    info <- c(
        info, 
        '##INFO=<ID=Check-Score,Number=1,Type=Float,Description="StemCNV Check-Score for CNV call">',     
        '##INFO=<ID=Precision,Number=1,Type=Float,Description="Estimated precision for this call">',
        '##INFO=<ID=HighImpact,Number=1,Type=String,Description="Overlapping high impact sites (StemCNV-check defined)">',
        '##INFO=<ID=Highlight,Number=1,Type=String,Description="Overlapping highlight sites (COSMIC genes)">',
        '##INFO=<ID=ROI,Number=1,Type=String,Description="Overlapping ROI sites (user defined)">',
        '##INFO=<ID=Gap_percent,Number=1,Type=Float,Description="Percent of segment which has a gap of probe coverage">',
        '##INFO=<ID=Genes,Number=1,Type=String,Description="Overlapping genes, sepearted by | character">'
    )
    }
    
    format <- c(
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Segment genotype">',
        '##FORMAT=<ID=CN,Number=1,Type=Integer,Description="Copy-number (estimated)">',
        '##FORMAT=<ID=TOOL,Number=1,Type=String,Description="Details for copy number calling tools">',
        '##FORMAT=<ID=LRR,Number=1,Type=Float,Description="Segment median Log R Ratio">'
        #FIXME (future): add clustered BAF
        # '##FORMAT=<ID=BAF,Number=1,Typeq=Foat,Description="Segment me(di)an B Allele Frequency">',
    )
    if (extra_annotation) {
    format <- c(format,
        '##INFO=<ID=REFCOV,Number=1,Type=Float,Description="Percentage of segment with matching call in reference sample">'
    )
    }
    # Filter is applied before extra annotation, so those calls are currently NOT in the VCF
    if (extra_annotation) { 
        density.perc.cutoff <- toolconfig$density.quantile.cutoff * 100 %>% round(1)
        gap.min.perc <- toolconfig$density.quantile.cutoff * 100 %>% round(1)
        filter <- c(
            '##FILTER=<ID=PASS,Description="All filters passed">',
            str_glue('##FILTER=<ID=high_probe_dens,Description="Probe density of segment is higher than {density.perc.cutoff}% of the array">'),
            str_glue('##FILTER=<ID=probe_gap,Description="Probe coverage of segment has considerbale gap (min. {gap.min.perc}%)">')
        )
    } else {
        filter.minsize <- toolconfig$filter.minsize
        filter.minprobes <- toolconfig$filter.minprobes
        filter.mindensity.Mb <- toolconfig$filter.mindensity
        filter <- c(
            '##FILTER=<ID=PASS,Description="All filters passed">',
            str_glue('##FILTER=<ID=Size,Description="CNV call <below min. size <{filter.minsize}bp">'),
            str_glue('##FILTER=<ID=n_probes,Description="CNV call from <{filter.minprobes} probes">'),
            str_glue('##FILTER=<ID=Density,Description="CNV call with <{filter.mindensity.Mb} probes/Mb">')
        )
    }
    
    out <- c()
    if (INFO) { out <- c(out, info) }
    if (FORMAT) { out <- c(out, format) }
    if (FILTER) { out <- c(out, filter) }
    return(out)
    
}

fix_header_lines <- function(header_lines, regex = NULL, contig_format = NULL) {
    
    if (!is.null(regex)) {
        header_lines <- header_lines %>% str_subset(regex)
    }
    
    if (!is.null(contig_format)) {
        header_tb <- tibble(line = header_lines) %>%
            mutate(
                is_contig = str_detect(line, '^##contig'),
                orig_order = row_number(),
                chrom_name = ifelse(is_contig, str_extract(line, '(?<=ID=)[^,]+'), NA_character_),
                fixed_chrom = NA_character_
            )
        # mapSeqlevels won't work proplerly if it is handed any NA values
        header_tb[header_tb$is_contig,]$fixed_chrom <- header_tb[header_tb$is_contig,]$chrom_name %>% 
            mapSeqlevels(contig_format)
        header_lines <- header_tb %>%
            filter(!is_contig | !is.na(fixed_chrom)) %>%
            mutate(
                fixed_chrom = factor(fixed_chrom, levels = genomeStyles('Homo_sapiens')[[contig_format]]),
                new_order = ifelse(is_contig, min(orig_order[is_contig]) - 1 + as.integer(fixed_chrom), orig_order),
                line = str_replace(line, '(?<=ID=)[^,]+', as.character(fixed_chrom))
            ) %>%
            arrange(new_order) %>%
            pull(line)
    }
    
    return(header_lines)
    
    
}

# Additions to Consider:
# Fuzzyness of CNV calls (based on closest probes?)
# Include CIPOS, CILEN
get_fix_section <- function(tb) {
    # Return a matrix/df with the following columns:
    # CHROM POS ID REF ALT QUAL FILTER INFO
    
    # expected minimal cols in tb:
    # seqnames, start, end, width, ID, CNV_caller, CNV_type, CN, sample_id, 
    #  n_initial_calls, initial_call_details, n_probes, n_uniq_probes, probe_density_Mb
    #  FILTER
    
    base_info_str <- 'END={end};SVLEN={width};SVCLAIM=D;N_PROBES={n_probes};N_UNIQ_PROBES={n_uniq_probes};PROBE_DENS={probe_density_Mb}'
    extra_info_str <- paste(
        base_info_str,
        'Check-Score={`Check-Score`};Precision={Precision_Estimate}',
        'HighImpact={high_impact_hits};Highlight={highlight_hits};ROI={ROI_hits}',
        'Gap_percent={percent_gap_coverage};Genes={overlapping_genes}',
        sep=';'
    )
    # Technically should check for all columns, but they come in a bundle
    use_info_str <- ifelse('Check-Score' %in% colnames(tb), extra_info_str, base_info_str)
    
    tb %>%
        mutate(
            # Need to use across + any_of to make this work if the columns aren't there
            # replace , by | for separator in INFO cols
            across(
                any_of(c("high_impact_hits", "highlight_hits", "ROI_hits", "overlapping_genes")),
                ~ str_replace_all(., ',', '|')
            ),
            # round numbers
            across(where(is.numeric), ~ round(., 3)),
            # convert NA or empty string to ".", all columns with possible NA to character
            across(
                any_of(c("Check-Score", "Precision_Estimate", "high_impact_hits", "highlight_hits",
                         "ROI_hits", "percent_gap_coverage", "overlapping_genes")),
                ~ ifelse(is.na(.) | . == "", '.', as.character(.))
            ),
                       
            # From VCF specs:
            # Note that for structural variant symbolic alleles, POS corresponds 
            # to the base immediately preceding the variant.
            POS = start - 1,
            REF = '.',
            #FIXME (future): CN>=4 should probably not be DUP but CNV accoring to newest VCF specs
            ALT = str_glue("<{CNV_type}>"),
            QUAL = '.',
            # END position of the longest variant described in this record. The END of each allele is defined as
            # <DEL>, <DUP>, <INV>, and <CNV> symbolic structural variant alleles:, POS + SVLEN.
            # in granges: width = end - start + 1; so this fits with corrected POS
            INFO = str_glue(use_info_str)
        ) %>%
        dplyr::rename(CHROM = seqnames) %>%
        select(CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO) %>%
        as.matrix()
}

get_gt_section <- function(tb, sample_sex){
    # Return a matrix with the following columns:
    # FORMAT, <sample>
    
    # min. expected cols in tb:
    # seqnames, start, end, width, ID, CNV_caller, CNV_type, CN, sample_id, 
    #  n_initial_calls, initial_call_details, n_probes, n_uniq_probes, probe_density_Mb
    #  FILTER, LRR
    
    # FORMAT keys: GT, CN, TOOL (str desc), LRR (median), [opt: REFCOV]
    # Future: BAF (no. of clusters?)
    
    sex_chroms <- get_sex_chroms(tb)
    use_cols <- c('GT', 'CN', 'TOOL', 'LRR')
    if ('reference_coverage' %in% colnames(tb)) {
        use_cols <- c(use_cols, 'REFCOV')
    }

    tb %>%
        rename_with(~str_replace(., 'reference_coverage', 'REFCOV')) %>%
        mutate(
            across(any_of('REFCOV'), ~ ifelse(is.na(.), '.', round(., 3) %>% as.character())),
            FORMAT = paste(use_cols, collapse = ':'),
            GT = case_when(
                # Male sex chroms have lower CN baseline
                # Male CN=1 will never be present (not a CNV, no LOh is called)
                sample_sex == 'm' & seqnames %in% sex_chroms & CN == 2 ~ '0/1',
                sample_sex == 'm' & seqnames %in% sex_chroms & CN > 2  ~ './.',
                # All other chromosomes
                # LOH
                CN == 2 ~ './.',
                # hom loss
                CN == 0 ~ '1/1',
                # single copy gain/loss
                CN %in% c(1, 3) ~ '0/1',
                # multi copy gain
                CN > 3 ~ './.',
            ),
            CN = as.character(CN),
            # initial_call_details; this will already carry further extra annotation if existing
            initial_call_details = ifelse(is.na(initial_call_details), '.', initial_call_details),
            TOOL = str_glue("caller={CNV_caller};n_initial_calls={n_initial_calls};initial_call_details={initial_call_details}"),
        ) %>%
        select(FORMAT, sample_id, all_of(use_cols)) %>%
        group_by(FORMAT, sample_id) %>%
        # paste columns together
        unite("value", all_of(use_cols), sep = ':', remove = TRUE) %>%
        pivot_wider(names_from = sample_id, values_from = value, values_fn = list) %>%
        unnest(cols = unique(tb$sample_id)) %>%
        as.matrix()
}


# Open question for writing combined CNV vcf:
# - how/which format to write back after CNV_processing
# - yes/no keep 'non-overlapped' calls
#    >> also ask Harald & Vale what they think here?

# Include ref sample in VCF?
# ##PEDIGREE=<ID=DerivedID,Original=OriginalID>
# use gvcf style otherwise