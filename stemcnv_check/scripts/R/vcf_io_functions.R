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
        as_granges(seqnames = CHROM, start = POS, width = 1)
    if (apply_filter) {
        vcf <- vcf %>% filter(FILTER == 'PASS')
    }
    return(vcf)
}

# parse_cnv_vcf

# CNV vcf writing
static_cnv_vcf_header <- function(toolconfig, INFO = TRUE, FORMAT = TRUE, FILTER= TRUE) { 
    filter.minsize <- toolconfig$filter.minsize
    filter.minprobes <- toolconfig$filter.minprobes
    filter.mindensity.Mb <- toolconfig$filter.mindensity
    
    info <- c(
        '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the longest variant described in this record">',
	    '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of structural variant">',
        '##INFO=<ID=SVCLAIM,Number=A,Type=String,Description="Claim made by the structural variant call. Valid values are D, J, DJ for abundance, adjacency and both respectively">',
        # Maybe use Number=. for these, if we want to allow multiple values after call overlaps
        '##INFO=<ID=N_PROBES,Number=1,Type=Integer,Description="Number of array probes in segment">',
        '##INFO=<ID=N_UNIQ_PROBES,Number=1,Type=Integer,Description="Number of unique array probe positions in segment">',
        '##INFO=<ID=PROBE_DENS,Number=1,Type=Float,Description="Density of Probes in segment (Probes / 10Mb)">'
    )
    format <- c(
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Segment genotype">',
        '##FORMAT=<ID=CN,Number=1,Type=Integer,Description="Copy-number (estimated)">',
        '##FORMAT=<ID=TOOL,Number=1,Type=String,Description="Details for copy number calling tools">',
        '##FORMAT=<ID=LRR,Number=1,Type=Float,Description="Segment median Log R Ratio">'
        # #TODO: not sure how to summarize this, maybe some simple clustering and # of clusters?
        # '##FORMAT=<ID=BAF,Number=1,Typeq=Foat,Description="Segment me(di)an B Allele Frequency">',
    )
    filter <- c(
        '##FILTER=<ID=PASS,Description="All filters passed">',
        str_glue('##FILTER=<ID=Size,Description="CNV call <below min. size <{filter.minsize}bp">'),
        str_glue('##FILTER=<ID=n_probes,Description="CNV call from <{filter.minprobes} probes">'),
        str_glue('##FILTER=<ID=Density,Description="CNV call with <{filter.mindensity.Mb} probes/Mb">')
    )
    
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
get_fix_section <- function(tb){
    # Return a matrix/df with the following columns:
    # CHROM POS ID REF ALT QUAL FILTER INFO
    
    # expected cols in tb:
    # seqnames, start, end, width, ID, CNV_caller, CNV_type, CN, sample_id, 
    #  n_initial_calls, initial_IDs_with_CN, n_probes, n_uniq_probes, probe_density_Mb
    #  FILTER
    
    tb %>%
        mutate(
            # round floats
            probe_density_Mb = round(probe_density_Mb, 2),
            # From VCF specs:
            # Note that for structural variant symbolic alleles, POS corresponds 
            # to the base immediately preceding the variant.
            POS = start - 1,
            REF = '.',
            #TODO: CN>=4 should probably not be DUP but CNV accoring to newest VCF specs
            ALT = str_glue("<{CNV_type}>"),
            QUAL = '.',
            # END position of the longest variant described in this record. The END of each allele is defined as
            # <DEL>, <DUP>, <INV>, and <CNV> symbolic structural variant alleles:, POS + SVLEN.
            # in granges: width = end - start + 1; so this fits with corrected POS
            INFO = str_glue(
                'END={end};SVLEN={width};SVCLAIM=D;N_PROBES={n_probes};N_UNIQ_PROBES={n_uniq_probes};PROBE_DENS={probe_density_Mb}'
            )
        ) %>%
        dplyr::rename(
            CHROM = seqnames,
        ) %>%
        select(CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO) %>%
        as.matrix()
}

get_gt_section <- function(tb, snp_vcf_gr){
    # Return a matrix/df with the following columns:
    # FORMAT, sample1, ...
    
    # expected cols in tb:
    # seqnames, start, end, width, ID, CNV_caller, CNV_type, CN, sample_id, 
    #  n_initial_calls, initial_call_details, n_probes, n_uniq_probes, probe_density_Mb
    #  FILTER
    
    # FORMAT keys: GT, CN, TOOL (str desc), LRR (median)
    # Future: BAF (no. of clusters?)
    
    median_lrr <- snp_vcf_gr %>%
        plyranges::select(LRR) %>%
        join_overlap_inner(as_granges(tb)) %>%
        as_tibble() %>%
        group_by(ID) %>%
        # round to 3 decimals
        summarise(LRR = median(LRR) %>% round(3))
    
    tb %>%
        select(sample_id, ID, CN, CNV_caller, n_initial_calls, initial_call_details) %>%
        mutate(
            FORMAT = 'GT:CN:TOOL:LRR',
            CN = as.character(CN),
            #TODO: X&Y on male will be different!
            GT = case_when(
                # LOH
                CN == 2 ~ './.',
                # hom loss
                CN == 0 ~ '1/1',
                # single copy gain/loss
                CN %in% c(1, 3) ~ '0/1',
                # multi copy gain
                CN > 3 ~ './.',
            ),
            # initial_call_details
            initial_call_details = ifelse(is.na(initial_call_details), '.', initial_call_details),
            TOOL = str_glue("caller={CNV_caller};n_initial_calls={n_initial_calls};initial_call_details={initial_call_details}"),
        ) %>%
        left_join(median_lrr, by = 'ID') %>%
        # Future TODO: some way to get (PennCNV) GC correction for LRR?
        # Future TODO: add solution for BAF
        select(FORMAT, sample_id, GT, CN, TOOL, LRR) %>%
        group_by(FORMAT, sample_id) %>%
        reframe(value = paste(GT, CN, TOOL, LRR, sep = ':')) %>%
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