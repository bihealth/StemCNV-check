suppressMessages(library(tidyverse))
suppressMessages(library(yaml))

get_sample_SNV_tb <- function(
    sample_gr, roi_tb, SNV_hotspot_table, gtf_file, ginfo_file, target_chrom_style, config
) {
    
    hotspot_genes <- c(
        roi_tb %>% filter(mapping == 'gene_name') %>% pull(hotspot),
        SNV_hotspot_table$gene_name
    ) %>% unique()
    gr_genes <- load_gtf_data(gtf_file, config, target_chrom_style, include_hotspot_genes = hotspot_genes)
    gr_info  <- load_genomeInfo(ginfo_file, config, target_chrom_style)
    
    #Note: this should ideally be parsed from vcf header
    #Note: can we get rsIDs from here also? (they are not consistently used as SNP probe IDs)
    meahri_ann_header <- c(
        'Allele', 'Annotation', 'Annotation_Impact', 'gene_name', 'gene_id', 'Feature_Type', "Feature_ID",
        'Transcript_BioType', 'Rank', 'HGVS.c', 'HGVS.p', 'cDNA.pos', 'CDS.pos', 'AA.pos',        
        # cDNA.pos / cDNA.length | CDS.pos / CDS.length | AA.pos / AA.length |
        'Distance', 'Strand', 'errors' # ERRORS / WARNINGS / INFO
    )
    sample_SNV_tb <- sample_gr %>%
        annotate_roi(roi_tb, gr_genes, gr_info, config) %>%
        as_tibble() %>%
        separate(ANN, sep = '\\|', into = meahri_ann_header) %>%
        dplyr::rename(
            GenCall_Score = IGC,
            Transcript_ID = Feature_ID
        ) %>%
        mutate(
            GT = ifelse(is.na(GT), './.', GT)
        )
    
    sample_SNV_tb
}


sample_GT_distances <- function(sample_SNP_gr, ref_SNP_gr, extra_snp_files, ref_SNP_vcf_file, target_chrom_style, use_filter=TRUE) {
    
    if (length(ref_SNP_gr) > 0) {
        # remove the ref_sample from extra_snp_files, is is already loaded
        extra_snp_files <- extra_snp_files[extra_snp_files != ref_SNP_vcf_file]
    }
    
    # Reading the SNP vcf files is the slowest thing here
    # > parallelized via furrr:future_map to speed it up, esp. for large number of files
    # FIXME (future): maybe don't parse the actual VCF here and use a faster table reader (readr / data.table) instead?
    # Or: The VariantAnnotation package might also have a faster vcf parser function
    furrr::future_map(
        extra_snp_files,
        \(x) {
            parse_snp_vcf(x, format_fields = c('GT'), apply_filter = use_filter) %>%
                fix_CHROM_format(target_chrom_style)
        },
        .options = furrr_options(seed = TRUE)
    ) %>%
        bind_ranges(sample_SNP_gr, ref_SNP_gr) %>%
        select(ID, sample_id, GT) %>%
        as_tibble() %>%
        mutate(GT = ifelse(GT == './.', NA, str_count(GT, '1'))) %>%
        pivot_wider(names_from = sample_id, values_from = GT, values_fill = NA) %>%
        filter(if_all(everything(), ~!is.na(.))) %>%
        dplyr::select(-seqnames, -start, -end, -strand, -width, -ID) %>%
        t() %>%
        dist(method = 'manhattan') %>%
        as.matrix() %>% 
        as.data.frame() %>%
        as_tibble(rownames = 'sample_distance_to')
}



get_SNV_table <- function(
    sample_SNV_tb,
    ref_SNP_gr,
    SNV_hotspot_table,
    config,
    defined_labels
) {
    
    subconfig <- config$settings$SNV_analysis
    
    ref_tb <- ref_SNP_gr %>%
        as_tibble() %>%
        dplyr::rename(
            ref_GenCall_Score = IGC,
            ref_GT = GT
        ) %>%
        group_by(seqnames, start, REF, ALT) %>%
        # Note: somehow slice_max is much slower than summarise/reframe ?
        # slice_max(ref_GenCall_Score, n=1, with_ties = FALSE) %>%
        summarise(
            # pick the GT per position with the highest GenCall_Score, ignore NA genotypes
            ref_GT = ifelse(
                n_distinct(ref_GT) == 1 | all(is.na(ref_GenCall_Score)),
                unique(ref_GT),
                unique(ref_GT[ref_GenCall_Score == max(ref_GenCall_Score)]) %>% na.omit()
            ),
            # max causes warnings if the input vector is zero-length
            ref_GenCall_Score = suppressWarnings(max(ref_GenCall_Score)),
            n_distinct_ref_GT = n_distinct(ref_GT),
        ) %>%
        mutate(ref_GT = ifelse(is.na(ref_GT), './.', ref_GT))
    
    if (any(ref_tb$n_distinct_ref_GT > 1)) {
        ref_tb <- ref_tb %>%
            group_by(seqnames, start, REF, ALT) %>% 
            summarise(
                ref_GT = paste(ref_GT, collapse = ' ; '),
                ref_GenCall_Score = max(ref_GenCall_Score),
            )
        warning(paste(
            'Multiple different GT calls with the same GenCall score for the same position in the reference sample.',
            'Collapsing to one GT call, this will not match the sample call:',
            ref_tb %>% filter(str_detect(ref_GT, ';')) %>%
                mutate(desc = str_glue('{seqnames}:{start}:{REF}>{ALT} - {ref_GT}')) %>% pull(desc) %>% paste(collapse = '\n')
        ))
    } else {
        ref_tb <- ref_tb %>% select(-n_distinct_ref_GT)
    }
    
    #Use whitespace as dummy regex that should match nothing (str_detect can not handle empty regex)
    protein_ablation_regex <- ifelse(
        is.null(subconfig$protein_ablation_annotations$Annotation_regex),
        ' ',
        subconfig$protein_ablation_annotations$Annotation_regex
    )
    protein_change_regex <- ifelse(
        is.null(subconfig$protein_change_annotations$Annotation_regex),
        ' ',
        subconfig$protein_change_annotations$Annotation_regex
    )
    selection_regex <- ifelse(
        is.null(subconfig$variant_selection$Annotation_regex),
        ' ',
        subconfig$variant_selection$Annotation_regex
    )
    snv_tb <- sample_SNV_tb %>%
        filter(!str_detect(FILTER, 'FEMALE-Y')) %>%
        # Only keep variants with:
        # - ROI overlap (if enabled)
        # - impact listed in any of the config sections
        # - annotation string matching any Annotation regex
        filter(
            (!is.na(ROI_hits) & subconfig$variant_selection$include_all_ROI_overlaps) |
                Annotation_Impact %in% subconfig$variant_selection$Impact | 
                str_detect(Annotation, selection_regex)                
        ) %>%
        select(
            seqnames, start, REF, ALT, ID, FILTER, GT, GenTrain_Score, GenCall_Score, 
            Annotation, Annotation_Impact, gene_name, Transcript_ID, Transcript_BioType, 
            HGVS.c, HGVS.p, ROI_hits
        ) %>%
        dplyr::rename(Impact = Annotation_Impact) %>%
        # Collapse information from probes with non-unique positions
        # could also do a majroity vote here, but highest score should be enough
        group_by(pick(1:4)) %>%
        slice_max(GenCall_Score, with_ties = FALSE) %>%
        # Consider: maybe also include non-variant sites not matching ref (i.e. changes 'back' to WT)
        # > these are more than likely mis-calls in the reference though
        filter(str_detect(GT, '1')) %>%
        ungroup() %>%
        # add/merge in reference GT
        left_join(ref_tb, by = c('seqnames', 'start', 'REF', 'ALT')) %>%
        # Assign SNV categories & labels (see also: label_name_definitions.yaml)
        mutate(
            SNV_category = case_when(
                !is.na(ROI_hits)                                               ~ 'ROI-overlap',
                # SNV_hotspot_table$hotspot does NOT contain the "p." part of the HGVS.p notation
                paste0(gene_name, '::', str_remove(HGVS.p, '^p.')) %in% 
                    SNV_hotspot_table$hotspot                                  ~ 'hotspot-match',
                gene_name %in% SNV_hotspot_table$gene_name                     ~ 'hotspot-gene',
                (str_detect(Annotation, protein_ablation_regex) | 
                    Impact %in% subconfig$protein_ablation_annotations$Impact) ~ 'protein-ablation',
                (str_detect(Annotation, protein_change_regex) | 
                    Impact %in% subconfig$protein_change_annotations$Impact)   ~ 'protein-changing',
                TRUE                                                           ~ 'other'
            ) %>%
                factor(levels = names(defined_labels$SNV_categories)),
            
            SNV_label = case_when(
                GT == ref_GT                                                   ~ 'reference genotype',
                SNV_category %in% c(subconfig$critical_SNV, subconfig$reportable_SNV) & 
                    (
                        GenCall_Score < subconfig$flag_GenCall_minimum |
                        ref_GenCall_Score < subconfig$flag_GenCall_minimum
                    )                                                          ~ 'unreliable critical/reportable',
                SNV_category %in% subconfig$critical_SNV                       ~ 'critical',
                SNV_category %in% subconfig$reportable_SNV                     ~ 'reportable',
                TRUE                                                           ~ 'de-novo SNV'
            ) %>%
                 factor(levels = names(defined_labels$SNV_labels)),
        ) %>%
        arrange(SNV_label, SNV_category, seqnames, start)

    # Ensure that only the defined labels are used
    stopifnot(
        all(snv_tb$SNV_label %in% names(defined_labels$SNV_labels)) & !any(is.na(snv_tb$SNV_label)),
        all(snv_tb$SNV_category %in% names(defined_labels$SNV_categories)) & !any(is.na(snv_tb$SNV_category))
    )
    
    snv_tb
}

get_SNV_hotspot_coverage <- function(
    sample_SNV_tb,
    SNV_hotspot_table
) {

    get_coverage_str <- function(pos) {
        pos <- pos[pos != '']
        if (length(pos) == 0) return('0% (0/NA)')
        total <- str_extract(pos, '(?<=/)[0-9]+$') %>% as.integer() %>% unique()
        stopifnot(length(total) == 1)
        n_covered <- str_extract(pos, '^[0-9]+') %>% 
            as.integer() %>% .[. > 0] %>% unique() %>% length()
        str_glue('{round(100*n_covered/total, 2)}% ({n_covered}/{total})') %>%
            as.character()
    }
    
    # retrun early and prevent error if no SNVs are present
    if (nrow(sample_SNV_tb) == 0) {
        return(tibble(
            gene_name = NA_character_,
            hotspots = NA_character_,
            Transcript_ID = NA_character_,
            Transcript_BioType = NA_character_,
            cDNA_covered = NA_character_,
            CDS_covered = NA_character_,
            AA_covered = NA_character_
        ))
    }    
    
    # This will only find genes based on existing annotations from the SNP vcf
    # therefore genes/hotspots not covered there will not be missing
    gene_coverage_table <- sample_SNV_tb %>%
        filter(!str_detect(FILTER, 'FEMALE-Y')) %>%
        filter(gene_name %in% SNV_hotspot_table$gene_name) %>%    
        group_by(gene_name) %>%
        group_modify(~ {
            spots <- SNV_hotspot_table %>% filter(gene_name == unique(.y$gene_name)) %>%
                filter(!is.na(HGVS.p)) %>% pull(HGVS.p)
            if (length(spots) == 0) {
                .x$hotspots <- 'none defined'
                return(.x)
            }
            covered <- spots[spots %in% .x$HGVS.p]
            not_covered <- spots[!spots %in% .x$HGVS.p]

            mutate(.x, 
                hotspots = paste0(
                    'covered (', length(covered), '/', length(spots), '): ', str_c(covered, collapse = ', '), '\n',
                    'missing (', length(not_covered), '/', length(spots), '): ', str_c(not_covered, collapse = ', ')
                )
            )
        }) %>%
        group_by(gene_name, hotspots, Transcript_ID, Transcript_BioType) %>%
        summarise(
            # Note: mehari  cDNA & CDS, and additionally UTR variants from CDS
            cDNA_covered = get_coverage_str(cDNA.pos[Distance==0]), #[!str_detect(HGVS.c, 'c\\.[0-9]+[+-][0-9]+[ACGT]>[ACGT]')]
            CDS_covered = get_coverage_str(CDS.pos[Distance==0]), #[str_detect(HGVS.c, 'c\\.[0-9]+[ACGT]>[ACGT]')]
            AA_covered = get_coverage_str(AA.pos),
        ) %>%
        group_by(gene_name, hotspots) %>%
        # select 1 transcript per gene
        arrange(
            gene_name,
            desc(str_detect(Transcript_BioType, 'ManeSelect')),
            desc(str_detect(Transcript_BioType, 'Coding')),
            desc(str_extract(AA_covered, '^[0-9.]+') %>% as.numeric()),
            desc(str_extract(CDS_covered, '(?<=/)[0-9]+(?=\\)$)') %>% as.numeric()),
        ) %>%
        slice(1)
    
    # Create summary for missing genes
    missing_genes <- SNV_hotspot_table %>%
        filter(!gene_name %in% gene_coverage_table$gene_name) %>%
        group_by(gene_name) %>%
        summarise(
            hotspots = ifelse(
                dplyr::n() == 1 & all(is.na(HGVS.p)), 
                'none defined',
                paste0(
                    'covered (0/', dplyr::n(), '): \n',
                    'missing (', dplyr::n(), '/', dplyr::n(), '): ', str_c(HGVS.p, collapse = ', ')
                )
            ),
            Transcript_ID = NA_character_,
            Transcript_BioType = NA_character_,
            cDNA_covered = get_coverage_str(''),
            CDS_covered = get_coverage_str(''),
            AA_covered = get_coverage_str(''),
        )


    bind_rows(
        gene_coverage_table,
        missing_genes
    ) %>%
        ungroup() %>%
        arrange(gene_name)
}


get_SNV_QC_table <- function(sample_id, sample_SNV_tb, ref_SNP_gr, SNV_table, use_filter) {
    
    # Get SNV/SNP summary values
    if (length(ref_SNP_gr) > 0) {
        format_size <- function(size) {
            oom_f <- case_when(
                size < 1e3 ~ 1,
                size < 1e6 ~ 1e3,
                size >= 1e6 ~ 1e6
            )
            string <- c('1' = 'bp', '1000' = 'kbp', '1e+06' = 'Mbp')
            paste0(round(size / oom_f), string[as.character(oom_f)])
        }
        
        snv_qc_tb <- sample_SNV_tb %>%
            left_join(
                ref_SNP_gr %>%
                    select(ID, FILTER, GT) %>%
                    as_tibble() %>%
                    dplyr::rename(ref_GT = GT, ref_FILTER = FILTER),
                by = c('seqnames', 'start', 'end', 'width', 'strand', 'ID')
            ) %>%
            mutate(
                numeric_GT = ifelse(GT == './.', NA, str_count(GT, '1')),
                numeric_GT_ref = ifelse(ref_GT == './.', NA, str_count(ref_GT, '1')),
                pw_filter_apply = !use_filter | (FILTER == 'PASS' & ref_FILTER == 'PASS')
            ) %>%
            filter(GT != ref_GT & !is.na(GT) & !is.na(ref_GT) & pw_filter_apply) %>%
            as_tibble() %>%
            group_by(sample_id, seqnames) %>%
            summarise(
                GT_distance = sum(
                    abs(numeric_GT - numeric_GT_ref),
                    na.rm = TRUE
                ),
                span = max(end) - min(start),
            ) %>%
            # Keep per chromsome info
            mutate(
                SNP_pairwise_distance_to_reference = sum(GT_distance),                
                chromosome = paste0(
                    case_when(
                        GT_distance == max(GT_distance) & span == min(span) ~ '[Most changes, Shortest span] ',
                        GT_distance == max(GT_distance) ~ '[Most changes] ',
                        span == min(span) ~ '[Shortest span] ',
                        TRUE ~ ''
                    ),                    
                    seqnames, ': ',
                    GT_distance, ' SNPs (over ',
                    format_size(span), ')'                    
                )
            ) %>%
            select(
                sample_id,
                SNP_pairwise_distance_to_reference,
                chromosome
            )
    } else {
        snv_qc_tb <- tibble(
            sample_id = sample_id,            
            SNP_pairwise_distance_to_reference = NA_real_,
            chromosome = NA_character_
            # chromosomes_with_changes = NA_character_,
            # chromsome_with_most_changes = NA_character_,
            # chromsome_with_shortest_span = NA_character_
        )
    }
    
    # number of reportable & critical SNVs
    snv_qc_tb$reportable_SNVs <- SNV_table %>%
        filter(SNV_label == 'reportable') %>%
        nrow()
    snv_qc_tb$critical_SNVs <- SNV_table %>%
        filter(SNV_label == 'critical') %>%
        nrow()
    
    # post filter SNPs
    n_filter_passed <- sum(sample_SNV_tb$FILTER == 'PASS')
    snv_qc_tb$SNPs_post_filter <- paste0(
        round(n_filter_passed/nrow(sample_SNV_tb) * 100, 2), '% (', n_filter_passed, ' SNPs)'
    )
    
    snv_qc_tb
}