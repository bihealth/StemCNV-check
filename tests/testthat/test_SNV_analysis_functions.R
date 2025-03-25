library(tidyverse)
library(plyranges)
library(testthat)
library(yaml)

source(test_path("../../src/stemcnv_check/scripts/R/helper_functions.R"))
source(test_path("../../src/stemcnv_check/scripts/R/hotspot_functions.R"))
source(test_path("../../src/stemcnv_check/scripts/R/vcf_io_functions.R"))
source(test_path("../../src/stemcnv_check/scripts/R/CNV_annotation_functions.R"))
source(test_path("../../src/stemcnv_check/scripts/R/snv_analysis_functions.R"))

# Functions to test:
# - [x] get_sample_SNV_tb
# - [x] get_SNV_table
# - [x] get_SNV_hotspot_coverage
# - [ ] get_SNV_QC_table
#      >> pw distance: 
# - [ ] sample_GT_distances
#      >> requires multiple (min. 3) vcf files, all need ID & GT set to usable values

config <- list(
    sample_table = test_path('../data/sample_table.xlsx'),
    snakedir = '',
    settings = list(
        SNV_analysis = list(
            snv_hotspot_table = test_path('../data/minimal-snv-hotspots.tsv'),
            flag_GenCall_minimum = 0.2,
            variant_selection = list(
              Impact = list('HIGH', 'MODERATE'),
              Annotation_regex = NULL,
              include_all_ROI_overlaps = TRUE
            ),
            critical_SNV = list('ROI-overlap', 'hotspot-match'),
            reportable_SNV = list('hotspot-gene', 'protein-ablation'),
            protein_ablation_annotations = list(
                Impact = list('HIGH'),
                Annotation_regex = NULL
            ),
            protein_change_annotations = list(
                Impact = list(),
                Annotation_regex = 'missense_variant|inframe'
            )
        ),
        CNV_processing = list(
            'gene_overlap' = list(
                'exclude_gene_type_regex' = c(),
                'include_only_these_gene_types' = c('lncRNA', 'miRNA', 'protein_coding'),
                'whitelist_hotspot_genes' = FALSE,
                'stemcell_hotspot_list' = test_path('../data/minimal-hotspots.tsv'),
                'cancer_gene_list' = test_path('../data/minimal-hotspots.tsv')
            ),
            'Check_score_values' = list('any_roi_hit' = 50)
        ),
        'vcf_output' = list(
            'chrom_style' = 'keep-original'
        )
    )
)

defined_labels <- yaml.load_file(test_path('../../src/stemcnv_check/control_files/label_name_definitions.yaml'))

SNV_hotspot_table <- load_hotspot_table(config, 'snv_hotspot')

empty_roi_tb <- tibble(
    list_name = character(),
    hotspot = character(),
    mapping = character(),
    call_type = character(),
    check_score = numeric(),
    description = character(),
    description_doi = character(),
    description_htmllinks = character(),
    order = numeric()
)

roi_tb <- tibble(
    list_name = 'ROI',
    hotspot = c('1:5000-5250', '1:8000-8250', 'DDX11L1'),
    mapping = c('position', 'position', 'gene_name'),
    call_type = 'any',
    check_score = 50,
    description = c('ROI_1: 1:5000-5250', 'ROI_2: 1:5000-5250', 'DDX11L1'),
    description_doi = NA_character_,
    description_htmllinks = c('ROI_1: 1:5000-5250', 'ROI_2: 1:5000-5250', 'DDX11L1'),
    order = 1:3,
)

# sample_SNV_tb from fixtures
sample_SNV_tb_roi <- sample_SNV_tb %>%
    mutate(ROI_hits = c(rep(NA, 20), '1:5000-5250', rep(NA, 9), '1:8000-8250', rep(NA, 80)))

#get_sample_SNV_tb(sample_gr, sample_id, SNV_hotspot_table, gtf_file, ginfo_file, target_chrom_style, config) 
test_that('get_sample_SNV_tb', {

    sample_id <- 'annot_sample'
    gtf_file <- test_path('../data/hg_minimal.gtf')
    ginfo_file <- test_path('../data/gr_info_minimal.tsv')
    target_chrom_style <- 'UCSC'

    sample_SNP_gr <- parse_snp_vcf(
        test_path('../data/minimal_probes_annotated.vcf'),
        info_fields = c('GenTrain_Score', 'ANN'),
        format_fields = c('GT', 'IGC'),
        apply_filter = FALSE
    ) 
    
    get_sample_SNV_tb(sample_SNP_gr, empty_roi_tb, SNV_hotspot_table, gtf_file, ginfo_file, target_chrom_style, config) %>% 
        expect_equal(sample_SNV_tb)
    
    # test with ROI(s)
    get_sample_SNV_tb(sample_SNP_gr, roi_tb, SNV_hotspot_table, gtf_file, ginfo_file, target_chrom_style, config) %>%
        expect_equal(sample_SNV_tb_roi)
})


#get_SNV_hotspot_coverage(sample_SNV_tb, SNV_hotspot_table) 
test_that('get_SNV_hotspot_coverage', {   
    # Listed hotspots: 'dummy_2.1', 'dummy_3.2', 'dummy_4.4', 'dummy_hotspot'
    # dummy_2.1, 4 variants: 1 intronic, 3 different sites\
    # dummy_3.2 & dummy_4.4: no variants in test data
    # dummy_hotspot, 2 variants: 1 frameshift, 1 missense (same AA)

    # Note: example file does NOT contain cDNA/CDS positions for intron variants, which mehari does until v0.30-something
    # > 1) cDNA/CDS positions on i.e. introns will be fixed in mehari
    # > 2) need to decide how to deal with variants that have more than one REF or ALT base >> unclear how to count them
    
    expected_tb <- tibble(
        gene_name = c('dummy_2.1', 'dummy_3.2', 'dummy_4.4', 'dummy_hotspot'),
        hotspots = c('covered (1/2): p.Pro200Glu\nmissing (1/2): p.Ala100Val', 'covered (0/1): \nmissing (1/1): p.Ala100Val', 'none defined', 'covered (1/1): p.Leu20Val\nmissing (0/1): '),
        Transcript_ID = c('ENST00000654321', NA, NA, 'ENST00000654456'),
        Transcript_BioType = c('Coding', NA, NA, 'Coding'),
        cDNA_covered = c('0.2% (3/1500)', '0% (0/NA)', '0% (0/NA)', '0.2% (2/1000)'),
        CDS_covered = c('0.3% (3/999)', '0% (0/NA)', '0% (0/NA)', '0.6% (2/333)'),
        AA_covered = c('0.9% (3/333)', '0% (0/NA)', '0% (0/NA)', '0.9% (1/111)'),
    )

    get_SNV_hotspot_coverage(sample_SNV_tb, SNV_hotspot_table) %>% 
        expect_equal(expected_tb)
})

# Available genes & variants
# dummy_1: transcript|ENST00000123456|Noncoding; gtf gene_type "unprocessed_pseudogene"
#   IMPACT: MODIFIER
#   Annotation: 4x non_coding_transcript_exon_variant, 1x downstream_gene_variant
# dummy_2.1: transcript|ENST00000654321|Coding; gtf gene_type "something_else"
#   IMPACT: 1x HIGH, 2x MODERATE, 1x LOW/MODIFIER
#   Annotation: 1x stop_gained, 2x missense_variant, 1x intron_variant
# !! no other genes from the gtf have probes, no other genes listed here are in the gtf !!
# dummy_hotspot: transcript|ENST00000654456|Coding; [no gtf]
#   IMPACT: 1x HIGH, 1x MODERATE
#   Annotation: 1x frameshift_variant, 1x missense_variant

#get_SNV_table(sample_SNV_tb, ref_SNP_gr, SNV_hotspot_table, config, defined_labels)
test_that('get_SNV_table', {
    outcols <- c(
        'seqnames', 'start', 'REF', 'ALT', 'ID', 'FILTER', 'GT', 'GenTrain_Score', 'GenCall_Score',
        'Annotation', 'Impact', 'gene_name', 'Transcript_ID', 'Transcript_BioType',	'HGVS.c', 'HGVS.p', 'ROI_hits'
    )
    test_config <- config
    empty_ref_gr <- GRanges(
        REF = character(),
        ALT = character(),
        IGC = integer(),
        GT = character()
    )
    
    # Test defaults without ref
    expected_tb <- sample_SNV_tb %>%
        dplyr::rename(Impact = Annotation_Impact) %>%
        select(all_of(outcols)) %>%
        filter(Impact %in% c('HIGH', 'MODERATE')) %>%
        mutate(
            ref_GT = NA, # function will output ligcal NA vector
            ref_GenCall_Score = NA_real_,
            SNV_category = c('hotspot-gene', 'hotspot-gene', 'hotspot-match', 'hotspot-gene', 'hotspot-match') %>%
                factor(levels = defined_labels$SNV_category_labels),
            SNV_label = c('reportable', 'reportable', 'critical', 'reportable', 'critical') %>%
                factor(levels = defined_labels$SNV_labels),
        ) %>%
        arrange(SNV_label, SNV_category)
    get_SNV_table(sample_SNV_tb, empty_ref_gr, SNV_hotspot_table, config, defined_labels) %>%
        expect_equal(expected_tb)
    # test defaults with ref
    ref_SNP_gr <- parse_snp_vcf(
        test_path('../data/minimal_probes.vcf'),
        format_fields = c('GT', 'IGC'),
        apply_filter = FALSE
    ) %>% fix_CHROM_format('UCSC')
    expected_tb <- expected_tb %>%
        mutate(
            ref_GT = '0/0',
            ref_GenCall_Score = '0.8'
        )
    get_SNV_table(sample_SNV_tb, ref_SNP_gr, SNV_hotspot_table, config, defined_labels) %>%
        expect_equal(expected_tb)
    
    # test with ROI
    expected_tb <- sample_SNV_tb_roi %>%
        dplyr::rename(Impact = Annotation_Impact) %>%
        select(all_of(outcols)) %>%
        filter(Impact %in% c('HIGH', 'MODERATE') | !is.na(ROI_hits)) %>%
        mutate(
            ref_GT = NA, # function will output logical NA vector
            ref_GenCall_Score = NA_real_,
            SNV_category = c('ROI-overlap', 'ROI-overlap', 'hotspot-gene', 'hotspot-gene', 'hotspot-match', 'hotspot-gene', 'hotspot-match') %>%
                factor(levels = defined_labels$SNV_category_labels),
            SNV_label = c('critical', 'critical', 'reportable', 'reportable', 'critical', 'reportable', 'critical') %>%
                factor(levels = defined_labels$SNV_labels),
        ) %>%
        arrange(SNV_label, SNV_category)
    get_SNV_table(sample_SNV_tb_roi, empty_ref_gr, SNV_hotspot_table, test_config, defined_labels) %>%
        expect_equal(expected_tb)
    
    # > changed variant_selection
    test_config$settings$SNV_analysis$variant_selection$Impact <- c('MODIFIER')
    test_config$settings$SNV_analysis$variant_selection$include_all_ROI_overlaps <- FALSE
    expected_tb <- sample_SNV_tb_roi %>%
        dplyr::rename(Impact = Annotation_Impact) %>%
        select(all_of(outcols)) %>%
        filter(Impact == 'MODIFIER') %>%
        mutate(
            ref_GT = NA, # function will output logical NA vector
            ref_GenCall_Score = NA_real_,
            SNV_category = c(rep('other', 5), 'ROI-overlap', 'hotspot-gene') %>%
                factor(levels = defined_labels$SNV_category_labels),
            SNV_label = c(rep('de-novo SNV', 5), 'critical', 'reportable') %>%
                factor(levels = defined_labels$SNV_labels),
        ) %>%
        arrange(SNV_label, SNV_category)
    get_SNV_table(sample_SNV_tb_roi, empty_ref_gr, SNV_hotspot_table, test_config, defined_labels) %>%
        expect_equal(expected_tb)
    
    # > changed protein_ablation_annotations
    test_config$settings$SNV_analysis$variant_selection$Impact <- c('HIGH', 'MODERATE')
    test_config$settings$SNV_analysis$variant_selection$Annotation_regex <- 'exon_variant|gene_variant'
    test_config$settings$SNV_analysis$protein_ablation_annotations$Annotation_regex <- 'exon_variant'
    expected_tb <- sample_SNV_tb %>%
        dplyr::rename(Impact = Annotation_Impact) %>%
        select(all_of(outcols)) %>%
        filter(Impact %in% c('HIGH', 'MODERATE') | gene_name == 'dummy_1') %>%
        mutate(
            ref_GT = NA, # function will output logical NA vector
            ref_GenCall_Score = NA_real_,
            SNV_category = c(
                'other', rep('protein-ablation', 4), 'other',
                'hotspot-gene', 'hotspot-gene', 'hotspot-match', 'hotspot-gene', 'hotspot-match'
            ) %>%
                factor(levels = defined_labels$SNV_category_labels),
            SNV_label = c(
                'de-novo SNV', rep('reportable', 4), 'de-novo SNV',
                'reportable', 'reportable', 'critical', 'reportable', 'critical'
            ) %>%
                factor(levels = defined_labels$SNV_labels),
        ) %>%
        arrange(SNV_label, SNV_category)
    get_SNV_table(sample_SNV_tb, empty_ref_gr, SNV_hotspot_table, test_config, defined_labels) %>% 
        expect_equal(expected_tb)
    
    # > changed critical_SNV <> reportable_SNV
    test_config <- config
    test_config$settings$SNV_analysis$critical_SNV <- c('ROI-overlap', 'hotspot-match', 'hotspot-gene')
    expected_tb <- sample_SNV_tb %>%
        dplyr::rename(Impact = Annotation_Impact) %>%
        select(all_of(outcols)) %>%
        filter(Impact %in% c('HIGH', 'MODERATE')) %>%
        mutate(
            ref_GT = NA, # function will output logical NA vector
            ref_GenCall_Score = NA_real_,
            SNV_category = c('hotspot-gene', 'hotspot-gene', 'hotspot-match', 'hotspot-gene', 'hotspot-match') %>%
                factor(levels = defined_labels$SNV_category_labels),
            SNV_label = rep('critical', 5) %>% factor(levels = defined_labels$SNV_labels),
        ) %>%
        arrange(SNV_label, SNV_category)
    get_SNV_table(sample_SNV_tb, empty_ref_gr, SNV_hotspot_table, test_config, defined_labels) %>%
        expect_equal(expected_tb)
})