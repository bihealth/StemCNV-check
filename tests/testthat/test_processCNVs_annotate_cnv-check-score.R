library(testthat)

library(tidyverse)
library(plyranges)

source(test_path('../../stemcnv_check/scripts/R/R_io_functions.R'))
source(test_path('../../stemcnv_check/scripts/R/processCNVs_annotate_check-score.R'))

# Functions:
# - finalize_gr_to_tb # Note: will be relaced once VEP is in place
# - annotate_cnv.check.score
# - annotate_precision.estimates

config <- list(
  'static_data' = list(
    'genome_gtf_file' = test_path('../data/hg_minimal.gtf'),
    'genomeInfo_file' = test_path('../data/gr_info_minimal.tsv')
  ),
  'settings' = list(
    'CNV_processing' = list(
      'gene_overlap' = list(
        'exclude_gene_type_regex' = c(),
        'include_only_these_gene_types' = c()
      )
    ),
    'chromosomes' = paste0('chr', c(1:22, 'X', 'Y'))
  )
)

base_tb <- tibble(
    seqnames = 'chr1',
    start = c(4000, 10000, 28000000, 28060000, 40000, 5000000, 3000),
    end   = c(5500, 14000, 28055000, 28065000, 50000, 7000000, 60000),
    sample_id = 'test_sample',
    CNV_type = c('gain', 'gain', 'gain', 'loss', 'loss', 'LOH', 'LOH'),
    ID = paste('combined', CNV_type, seqnames, start, end, sep='_'),
    CNV_caller = list(c('toolA','toolB'), 'toolA', c('toolA','toolB'), c('toolA','toolA','toolB'), 'toolA', 'toolB', 'toolB'),
    #n_premerged_calls = list(c(2,1), 1, c(2,1), c(2,2,1),1,1,1),
    n_snp_probes = list(c(10,5), 100, c(100,50), c(100,100,50),100,50,50),
    n_uniq_probe_positions = c(15, 100, 150, 100, 100, 50, 50),
    copynumber = list(c('3','3'), '3', c('3','4'), c('1','1','1'), '1', '2', '2'),
    caller_confidence = list(c(1,1), 1, c(1,1), c(1,1,1), 1,1,1),
    #overlap_merged_call = NA_real_,
    caller_merging_coverage = c('toolA-100,toolB-100', NA, 'toolA-100,toolB-100', 'toolA-80,toolB-100', NA, NA, NA),
    caller_merging_state = c('combined', NA, rep('combined', 2), 'no-overlap', 'no-overlap', 'no-overlap'),
    reference_overlap = c(T, T, F, F, F, F, T),
    reference_coverage = c(100, 85.01, NA_real_, NA_real_, NA_real_, NA_real_, 60) %>% as.list(),
    reference_caller = list(c('toolA','toolB'), 'faketool', NA_character_, NA_character_,NA_character_,NA_character_, 'toolA'),
    high_impact_hits = c(NA, NA, 'dummyC', NA, '1p36,chr1:40000-50000', NA, NA),
    highlight_hits = c(NA, 'DDX11L1', NA, NA, NA, NA, 'DDX11L1,dummyB'),
    ROI_hits = c('fake-ROI', NA, NA, 'dummyC', NA, NA, NA),    
    percent_gap_coverage = c(0, 2000/4001, 25000/55001, 1000/5001, 2000/10001, 1e6/(2e6+1), 0),
    probe_coverage_gap = c(FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE),
    high_probe_density = c(NA, NA, TRUE, TRUE, NA, FALSE, FALSE)
)

expected_gr_to_tb <- base_tb %>%
  mutate(
    seqnames = factor(seqnames, levels = get_chromosome_set()),
    length = c( 1501, 4001, 55001, 5001, 10001, 2000001, 57001),
    n_premerged_calls = list(NULL),
    overlap_merged_call = NA_real_,
    n_genes = c(1, 1, 1, 1, 1, 0, 3),
    overlapping_genes = c('dummyA', 'DDX11L1', "dummyC", "dummyC", 'dummyB', NA, 'dummyA,DDX11L1,dummyB'),
    `Check-Score` = NA_real_,
    Precision_Estimate = NA_real_
  ) %>%
  # any_of vs one_of here, since we DON'T want to test yet 
  select(any_of(colnames(expected_final_tb)))

test_that("test gr_to_final_tb", {
    # Test scenarios:
    # - missing columns (list & normal)
    # - CNV has no gene overlap
    # - CNV has (partial) gene overlap
    # - CNV has multiple genes overlapping
    gr_genes <- load_gtf_data(config)
    expect_equal(finalise_gr_to_tb(as_granges(base_tb), gr_genes), expected_gr_to_tb)
    # extra tests:
    # - column should be list but is not
    expect_equal(finalise_gr_to_tb(base_tb %>%
                                     mutate(copynumber = c(3,3,3,1,1,2,2)) %>%
                                     as_granges(), 
                                   gr_genes), 
                 expected_gr_to_tb %>%
                   mutate(copynumber = c(3,3,3,1,1,2,2) %>% as.list())
    )
})


test_that("Annotate CNV check scores", {
  CNV_size_score <- function(len) {   1/3 * log(len) * log(len) - 15 }
  LOH_size_score <- function(len) { 0.275 * log(len) * log(len) - 15 }
  config$settings$CNV_processing$Check_score_values <- list(
      'roi_hit_base' = 50,
      'highimpact_base' = 20,
      'highlight_base' = 0,
      'per_gene_roi' = 10,
      'per_gene_highimpact' = 10,
      'per_gene_highlight' = 5,
      'per_gene_any' = 0.2
  )
  
  hi_gr <- tibble(
    seqnames = 'chr1',
    start = c(28050000, 40000, 0),
    end = c(28070000, 50000, 7200000),
    strand = c('+', '*', '*'),
    list_name = 'test-list',
    hotspot = c('dummyC', 'chr1:40000-50000', '1p36'),
    mapping = c('gene_name', 'position', 'gband'),
    call_type = c('gain', 'any', 'loss'),
    check_score = c(15, NA, NA),
    source = c('dummy', 'dummy', 'dummy'),
    comment = NA
  ) %>% as_granges()
  
  hl_gr <- tibble(
    seqnames = 'chr1',
    start = 11873,
    end = 14409,
    strand = '+',
    list_name = 'test-list',
    hotspot = 'DDX11L1', 
    mapping = 'gene_name', 
    call_type = 'any', 
    check_score = NA_real_,
    source = 'dummy',
    comment = NA
  ) %>% as_granges()
  
  expected_tb <- expected_gr_to_tb %>%
    mutate(
      # Test scenarios:
      `Check-Score` = c(
        # ROI (base + 1 hit) + 1 other gene 
        CNV_size_score(1501) + 50 + 10 + 0.2, 
        # HL hit (base is 0)
        CNV_size_score(4001) + 5, 
        # HI hit, 1 gene with custom score 15
        CNV_size_score(55001) + 20 + 15, 
        # ROI hit (base + 1), ROI *is* the one overlapping gene, so no extra
        CNV_size_score(5001) + 50 + 10,
        # HI hit (base  + 2 genes) , + 1 other gene
        CNV_size_score(10001) + 20 + 10+10 + 0.2,
        # no genes
        LOH_size_score(2000001) + 0,
        # HL hit (1 gene) + 2 other genes
        LOH_size_score(57001) + 5 + 0.4
      )
    )
  expect_equal(annotate_cnv.check.score(expected_gr_to_tb, hi_gr, hl_gr, config$settings$CNV_processing$Check_score_values), expected_tb)
  # Extra test:
  # - HL base score not 0
  config$settings$CNV_processing$Check_score_values$highlight_base <- 2.7
  expected_tb[c(2,7), 'Check-Score'] <- expected_tb[c(2,7), 'Check-Score'] + 2.7
  expect_equal(annotate_cnv.check.score(expected_gr_to_tb, hi_gr, hl_gr, config$settings$CNV_processing$Check_score_values), expected_tb)
})



#C)
# - CBS only
# - PennCNV only
# - CBS + PennCNV
# - different sizes