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

parse_snp_vcf <- function(
    vcf, info_fields = FALSE, format_fields = c("LRR", "BAF"), apply_filter = TRUE
) {
    if (typeof(vcf) == 'character') {
        vcf <- read.vcfR(vcf, verbose = F)
    }
    # VCF POS should be 1-based,
    # Granges are also 1-based
    # AND are (fully) open [= start & end are included]
    vcf <- vcf %>%
        vcfR_to_tibble(info_fields = info_fields, format_fields = format_fields) %>%
        as_granges(seqnames = CHROM, start = POS, width = 1)
    if (apply_filter) {
        vcf <- vcf %>% filter(FILTER == 'PASS')
    }
    return(vcf)
}

# parse_cnv_vcf
parse_cnv_vcf <- function(
    vcf, info_fields = NULL, format_fields = NULL, apply_filter = FALSE
) {
    if (typeof(vcf) == 'character') {
        vcf <- read.vcfR(vcf, verbose = F)
    }
    # vcfR2tidy does not work on empty vcf objects
    # FIXME: what about the additional annotation fields?
    if (nrow(vcf@fix) == 0) {
        vcf <- tibble(
            seqnames = character(),
            start = integer(),
            ID = character(),
            FILTER = character(),
            end = integer(),
            width = integer(),
            n_probes = integer(),
            n_uniq_probes = integer(),
            probe_density_Mb = double(),            
            sample_id = character(),
            # GT = character(),
            CN = integer(),
            LRR = double(),
            CNV_type = character(),
            CNV_caller = character(),
            n_initial_calls = integer(),
            initial_call_details = character(),
        ) %>%
            as_granges()
    } else {
        # VCF POS are 1-based, *but* for SV/CNV entries describe the base _before_ the variant
        # Granges are also 1-based and are (fully) open [= start & end should be included]
        vcf <- vcf %>%
            vcfR_to_tibble(info_fields = info_fields, format_fields = format_fields) %>%
            separate(TOOL, into = c('CNV_caller', 'n_initial_calls', 'initial_call_details'), sep = ';') %>%
            mutate(
                POS = POS + 1,
                ALT = str_remove(ALT, '<') %>% str_remove('>') %>%
                    str_replace('DUP', 'gain') %>%
                    str_replace('DEL', 'loss') %>%
                    str_replace('CNV:LOH', 'LOH'),
                across(c(CNV_caller, n_initial_calls, initial_call_details), ~ str_extract(., '(?<==)[^;]+')),
                # CNV_caller = str_extract(TOOL, '(?<=caller=)[^;]+'),
                # n_initial_calls = str_extract(TOOL, '(?<=n_initial_calls=)[^;]+'),
                # initial_call_details = str_extract(TOOL, '(?<=initial_call_details=)[^;]+'),
                FILTER = ifelse(FILTER %in% c('', 'PASS'), NA_character_, FILTER),
                across(
                    any_of('PREC_DESC'),
                    ~ str_replace_all(., '=', ': ') %>% str_replace_all(';', '; ')
                ),
                # FIXME: so far all gene & hotspots are NOT back_converted to comma-separated
                #  maybe do this and adapt report functions accordingly
                # across(any_of('overlapping_genes'), ~ str_replace_all(., '\\|', ', ')),
                across(
                    any_of(c('n_initial_calls', 'REFCOV', 'PREC_EST')),
                    ~ ifelse(. == '.', NA_real_, as.numeric(.))
                ),
                across(where(is.character), ~ ifelse(. == '.', NA_character_, .)),
            ) %>%
            rename_with(~str_to_lower(.), contains('PROBE')) %>%
            rename_with(~str_replace(., 'REFCOV', 'reference_coverage'), contains('REFCOV')) %>%
            rename_with(~str_replace(., 'ROI', 'ROI_hits'), contains('ROI')) %>%
            rename_with(
                ~str_replace(., 'PREC_DESC', 'precision_estimate_description') %>%
                    str_replace('PREC_EST', 'precision_estimate'),
                contains('PREC_')
            ) %>%
            dplyr::rename(CNV_type = ALT, probe_density_Mb = probe_dens) %>%
            select(-GT, -REF, -QUAL, -SVCLAIM) %>%
            as_granges(seqnames = CHROM, start = POS, end = END, width = SVLEN)
    }
        
    if (apply_filter) {
        vcf <- vcf %>% filter(is.na(FILTER))
    }
    return(vcf)
}


# CNV vcf writing
static_cnv_vcf_header <- function(
    toolconfig, extra_annotation = FALSE, INFO = TRUE, FORMAT = TRUE, FILTER = TRUE, ALT = FALSE, fullconfig = NULL
) {

    alt <- c('##ALT=<ID=CNV:LOH,Description="Loss of heterozygosity, same as run of homozygosity">')
    
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
        stopifnot(!is.null(fullconfig))
        info <- c(
            info,
            '##INFO=<ID=Check_Score,Number=1,Type=Float,Description="StemCNV Check_Score for CNV call">',
            paste0(
                '##INFO=<ID=Call_label,Number=1,Type=String,Description="Evaluation of CNV, based on reference ',
                'overlap, Check-Score and Filters; one of: ',
                paste(names(fullconfig$evaluation_settings$CNV_call_labels), collapse = ', '),
                ')">'
            ),
            '##INFO=<ID=stemcell_hotspot,Number=1,Type=String,Description="Overlapping stemcell hotspot sites (StemCNV-check defined)">',
            '##INFO=<ID=dosage_sensitive_gene,Number=1,Type=String,Description="Overlapping dosage sensitive genes (Collins et al. 2022)">',
            '##INFO=<ID=cancer_gene,Number=1,Type=String,Description="Overlapping cancer genes (Intogen cancer drivers)">',
            '##INFO=<ID=Gap_percent,Number=1,Type=Float,Description="Percent of segment which has a gap of probe coverage">',
            '##INFO=<ID=overlapping_genes,Number=1,Type=String,Description="Overlapping gene names">'
        )
    }
    
    format <- c(
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Segment genotype">',
        '##FORMAT=<ID=CN,Number=1,Type=Integer,Description="Copy-number (estimated)">',
        '##FORMAT=<ID=LRR,Number=1,Type=Float,Description="Segment median Log R Ratio">',
        '##FORMAT=<ID=TOOL,Number=1,Type=String,Description="Details for copy number calling tools">'
        #FIXME (future): add clustered BAF
        # '##FORMAT=<ID=BAF,Number=1,Typeq=Foat,Description="Segment me(di)an B Allele Frequency">',
    )
    if (extra_annotation) {
        min.ref.ov <- toolconfig$min.reciprocal.coverage.with.ref * 100 %>% round(1)
        format <- c(
            format,
            '##FORMAT=<ID=ROI,Number=1,Type=String,Description="Overlapping ROI sites">',
            '##FORMAT=<ID=PREC_EST,Number=1,Type=Float,Description="Estimated precision for this call">',
            '##FORMAT=<ID=PREC_DESC,Number=1,Type=String,Description="Data basis for precision esimate">',
            paste0(
                '##FORMAT=<ID=REFCOV,Number=1,Type=Float,Description="Percentage of segment with matching call ',
                str_glue(' in reference sample (min {min.ref.ov}% reciprocal overlap)">')
            )
        )
    }
    # See also label_name_definitions.yaml for all filters defined in the workflow 
    filter.minsize <- toolconfig$filter.minsize
    filter.minprobes <- toolconfig$filter.minprobes
    filter.mindensity.Mb <- toolconfig$filter.mindensity
    filter <- c(
        '##FILTER=<ID=PASS,Description="All filters passed">',
        str_glue('##FILTER=<ID=min_size,Description="CNV call <below min. size <{filter.minsize}bp">'),
        str_glue('##FILTER=<ID=min_probes,Description="CNV call from <{filter.minprobes} probes">'),
        str_glue('##FILTER=<ID=min_density,Description="CNV call with <{filter.mindensity.Mb} probes/Mb">')
    )
    if (extra_annotation) {
        density.perc.cutoff <- toolconfig$density.quantile.cutoff * 100 %>% round(1)
        gap.min.perc <- toolconfig$density.quantile.cutoff * 100 %>% round(1)
        filter <- c(
            '##FILTER=<ID=PASS,Description="All filters passed">',
            str_glue('##FILTER=<ID=high_probe_dens,Description="Probe density of segment is higher than {density.perc.cutoff}% of the array">'),
            str_glue('##FILTER=<ID=probe_gap,Description="Probe coverage of segment has considerbale gap (min. {gap.min.perc}%)">')
        )
    }
    
    out <- c()
    if (INFO) { out <- c(out, info) }
    if (FORMAT) { out <- c(out, format) }
    if (FILTER) { out <- c(out, filter) }
    if (ALT) { out <- c(out, alt) }
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
            group_by(is_contig) %>%
            mutate(
                fixed_chrom = factor(fixed_chrom, levels = genomeStyles('Homo_sapiens')[[contig_format]]),
                orig_order = ifelse(is_contig, rep(min(orig_order), dplyr::n()), orig_order),
                line = ifelse(is_contig, str_replace(line, '(?<=ID=)[^,]+', as.character(fixed_chrom)), line)
            ) %>%
            arrange(orig_order, fixed_chrom) %>%
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
        'Check_Score={Check_Score};Call_label={Call_label}',
        'stemcell_hotspot={stemcell_hotspot};dosage_sensitive_gene={dosage_sensitive_gene};cancer_gene={cancer_gene}',
        'Gap_percent={Gap_percent};overlapping_genes={overlapping_genes}',
        sep=';'
    )
    # Technically should check for all columns, but they come in a bundle
    use_info_str <- ifelse('Check_Score' %in% colnames(tb), extra_info_str, base_info_str)
    
    tb %>%
        arrange(seqnames, start) %>%
        mutate(
            # Need to use across + any_of to make this work if the columns aren't there
            # replace , by | for separator in INFO cols
            across(
                any_of(c("stemcell_hotspot", "dosage_sensitive_gene", "cancer_gene", "overlapping_genes")),
                ~ str_replace_all(., ',', '|')
            ),
            # round numbers
            across(where(is.numeric), ~ round(., 3)),
            # convert NA or empty string to ".", all columns with possible NA to character
            across(
                any_of(c(
                    "Check_Score", "Call_label", "stemcell_hotspot", "cancer_gene", "dosage_sensitive_gene",
                    "Gap_percent", "overlapping_genes"
                )),
                ~ ifelse(is.na(.) | . == "", '.', as.character(.))
            ),
            # Fix CNV_types used internally to match VCF specs
            # FIXME (future): CN>=4 should probably not be DUP but CNV accoring to newest VCF specs
            CNV_type = str_replace(CNV_type, 'gain', 'DUP') %>%
                str_replace('loss', 'DEL') %>%
                str_replace('LOH', 'CNV:LOH'),
            # From VCF specs:
            # Note that for structural variant symbolic alleles, POS corresponds 
            # to the base immediately preceding the variant.
            POS = as.character(start - 1),
            REF = '.',
            ALT = str_glue("<{CNV_type}>"),
            QUAL = '.',
            FILTER = ifelse(is.na(FILTER), 'PASS', FILTER),
            # END position of the longest variant described in this record. The END of each allele is defined as
            # <DEL>, <DUP>, <INV>, and <CNV> symbolic structural variant alleles:, POS + SVLEN.
            # in granges: width = end - start + 1; so this fits with corrected POS
            INFO = str_glue(use_info_str)
        ) %>%
        dplyr::rename(CHROM = seqnames) %>%
        select(CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO) %>%
        as.matrix()
}

get_gt_section <- function(tb, sample_id, sample_sex, target_style) {
    # Return a matrix with the following columns:
    # FORMAT, <sample>
    
    # sample_id is needed in case tb is empty (and/or doesn't contain any entries in sample_id column)
    if (nrow(tb) == 0) {
        out <- tibble(FORMAT = character())
        out[[sample_id]] <- character()
        return(out %>% as.matrix())
    }
    
    # min. expected cols in tb:
    # seqnames, start, end, width, ID, CNV_caller, CNV_type, CN, sample_id, 
    #  n_initial_calls, initial_call_details, n_probes, n_uniq_probes, probe_density_Mb
    #  FILTER, LRR
    # Additional columns used: reference_coverage, precision_estimate, precision_estimate_description
    
    # FORMAT keys: GT, CN, TOOL (str desc), LRR (median), [opt: ROI, PREC_EST, PREC_DESC, REFCOV]
    # Future: BAF (no. of clusters?)
    sex_chroms <- get_sex_chroms(target_style)
    use_cols <- c('GT', 'CN', 'LRR', 'TOOL')
    # Add extra columns from CNV processing, they should always be present in a set
    if ('reference_coverage' %in% colnames(tb)) {
        use_cols <- c(use_cols, 'ROI', 'PREC_EST', 'PREC_DESC', 'REFCOV')
    }

    tb %>%
        arrange(seqnames, start) %>%
        rename_with(
            ~str_replace(., 'reference_coverage', 'REFCOV') %>%
                str_replace('precision_estimate_description', 'PREC_DESC') %>%
                str_replace('precision_estimate', 'PREC_EST') %>%
                str_replace('ROI_hits', 'ROI')
        ) %>%
        mutate(
            across(any_of(c('PREC_EST', 'REFCOV')), ~ ifelse(is.na(.), '.', round(., 3) %>% as.character())),
            across(
                any_of(c('PREC_DESC', 'ROI')),
                ~ ifelse(is.na(.), '.', str_replace_all(., ': ', '=') %>% str_replace_all('; ', ';'))
            ),
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

# Include ref sample in VCF?
# ##PEDIGREE=<ID=DerivedID,Original=OriginalID>
# use gvcf style otherwise


# Also directly write out a cnv vcf
write_cnv_vcf <- function(cnv_tb, out_vcf, sample_sex, tool_name, fullconfig, snp_vcf_meta, command_desc, target_style) {

    stopifnot(str_ends(out_vcf, '\\.gz'))
    
    if (tool_name == 'combined_calls') {
        tool_config <- fullconfig$settings$CNV_processing$call_processing
    } else {
        tool_config <- fullconfig$settings[[tool_name]]
    }
        
    filtersettings <- tool_config$probe_filter_settings
    if (filtersettings == '_default_') {
        filtersettings <- fullconfig$settings$default_probe_filter_set
    }    
    # CBS does not have LOH calls, only combined calls need header info for extended annotation
    include_LOH_calls <- tool_name != 'CBS'
    extra_annotation <- tool_name == 'combined_calls'
    
    header <- c(
        # contains fileformat, so needs to go first
        fix_header_lines(snp_vcf_meta, 'fileformat|contig|BPM=|EGT=|CSV=', target_style),
        static_cnv_vcf_header(tool_config, extra_annotation = extra_annotation, ALT = include_LOH_calls, fullconfig = fullconfig),
        ifelse(str_starts(command_desc, '##'), command_desc, paste0('##', command_desc))
    )

    fix <- get_fix_section(cnv_tb)
    gt <- get_gt_section(cnv_tb, sample_id, sample_sex, target_style)
    # write.vcf does not work on empty vcfR objects
    if (nrow(cnv_tb) == 0) {
        vcf_file <- out_vcf %>% str_replace('.gz$', '')
        cat(header, file = vcf_file, sep = '\n')
        cat(
            paste0(
                '#', paste(c(colnames(fix), colnames(gt)), collapse = '\t')
            ),
            file = vcf_file, sep = '\n', append = T
        )
        R.utils::gzip(vcf_file)
    } else {
        cnv_vcf <- new(
            "vcfR",
            meta = header,
            fix = fix,
            gt = gt
        )
        write.vcf(cnv_vcf, out_vcf)
    }
}