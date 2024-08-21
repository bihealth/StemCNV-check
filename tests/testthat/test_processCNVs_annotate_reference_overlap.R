library(testthat)

library(tidyverse)
library(plyranges)

source(test_path('../../stemcnv_check/scripts/R/helper_functions.R'))
source(test_path('../../stemcnv_check/scripts/R/processCNVs_annotate_reference_overlap.R'))

# Example calls for testing

# Call with ref (100% match)
# Call with ref (85% match)
# Call with ref (wrong CNV type)
# Call with ref (ref <80% OV)
# Call with ref (call <80% OV)
# Call without ref
# Ref call without match

sample_cnvs <- tibble(   
  seqnames = 'chr1',
  start = c(4000, 10000, 12000, 15e4, 28e4, 35e4) %>% as.integer(),
  end   = c(5500, 11000, 14000, 20e4, 30e4, 40e4) %>% as.integer(),
  sample_id = 'test_sample',
  CNV_type = c('DUP', 'DUP', 'DEL', 'DEL', 'DEL', 'DEL'),
  ID = paste('combined', CNV_type, seqnames, start, end, sep='_'),
  CNV_caller = list(c('toolA','toolB'), 'toolA', c('toolA','toolB'), c('toolA','toolA','toolB'), 'toolA', 'toolB'),
  # n_premerged_calls = list(c(2,1), 1, c(2,1), c(2,2,1),1,1),
  n_probes = c(15, 5, 15, 25, 5, 5),
  CN = c(3, 3, 0, 1, 1, 1), 
      #list(c('3','3'), '3', c('1','0'), c('1','1','1'), '1', '1'),
  # caller_confidence = list(c(1,1), 1, c(1,1), c(1,1,1), 1,1),
  overlap_merged_call = NA_real_,
  caller_merging_coverage = c('toolA-100,toolB-100', NA, 'toolA-100,toolB-100', 'toolA-80,toolB-100', NA, NA),
  caller_merging_state = c('combined', NA, rep('combined', 2), 'no-overlap', 'no-overlap')
)

ref_cnvs <- tibble(   
  seqnames = 'chr1',
  start = c(4000, 10000, 12000, 18e4, 25e4, 45e4) %>% as.integer(),
  end   = c(5500, 10850, 14000, 20e4, 30e4, 50e4) %>% as.integer(),
  sample_id = 'test_sample',
  CNV_type = c('DUP', 'DUP', 'DUP', 'DEL', 'DEL', 'DEL'),
  ID = paste('combined', CNV_type, seqnames, start, end, sep='_'),
  CNV_caller = list(c('toolA','toolB'), 'faketool', c('toolA','toolB'), c('toolA','toolA','toolB'), 'toolA', 'toolB'),
  # n_premerged_calls = list(c(2,1), 5, c(2,1), c(2,2,1),1,1),
  n_probes = c(15, 20, 15, 25, 5, 5),
  CN = c(3, 4, 3, 1, 1, 1),
  # caller_confidence = list(c(1,1), 5, c(1,1), c(1,1,1), 1,1),
  overlap_merged_call = NA_real_,
  caller_merging_coverage = c('toolA-100,toolB-100', 'faketool', 'toolA-100,toolB-100', 'toolA-80,toolB-100', NA, NA),
  caller_merging_state = c(rep('combined', 4), 'no-overlap', 'no-overlap')
) %>% 
  bind_rows(get_expected_final_tb('UCSC')) %>%
  dplyr::select(-width, -reference_overlap, -reference_coverage, 
                -reference_caller, -n_genes, -overlapping_genes) %>%
  as_granges()

test_that("Annotate CNVs with ref", {
  #Note: this test will fail if `sample_cnvs` is converted to a granges object first and then
  # (changed back &) mutated, because somehow that makes the list columns appear as class 'AsIs'
  min.reciprocal.coverage.with.ref <- 80
  
  expected_gr <- sample_cnvs %>%
    mutate(
        reference_overlap = c(T, T, F, F, F, F),
        reference_coverage = c(100, 85.01, NA_real_, NA_real_, NA_real_, NA_real_),
        reference_caller = c('toolA;toolB', 'faketool', NA_character_, NA_character_,NA_character_,NA_character_)
    ) %>%
    as_granges()
  
  annotate_reference_overlap(
      as_granges(sample_cnvs),
      ref_cnvs,
      min.reciprocal.coverage.with.ref
  ) %>% 
    expect_equal(expected_gr) 
} )

test_that("Annotate CNVs empty ref", {
  min.reciprocal.coverage.with.ref <- 80
  
  expected_gr <- sample_cnvs %>%
    # This is somehow needed to successfully run the test via Rscript
    as_granges() %>% as_tibble() %>%
    mutate(
        reference_overlap = F, 
        reference_coverage = rep(NA_real_, 6),
        reference_caller = rep(NA_character_, 6)
    ) %>%
    as_granges()
  
  annotate_reference_overlap(
      as_granges(sample_cnvs), 
      ref_cnvs %>% filter(sample_id == 'non_existent_sample'),
      min.reciprocal.coverage.with.ref
  ) %>%
    expect_equal(expected_gr)        
} )

test_that("Annotate CNVs with empty input", {
  min.reciprocal.coverage.with.ref <- 80
  
  expected_gr <- sample_cnvs %>%
    # This is somehow needed to successfully run the test via Rscript
    as_granges() %>% as_tibble() %>%
    mutate(
        reference_overlap = F, 
        reference_coverage = rep(NA_real_, 6),
        reference_caller = rep(NA_character_, 6)
    ) %>%
    as_granges() %>%
    filter(sample_id == 'non_existent_sample')
  
  annotate_reference_overlap(
      sample_cnvs %>% as_granges() %>% filter(sample_id == 'non_existent_sample'),
      ref_cnvs, 
      min.reciprocal.coverage.with.ref
  ) %>%
    expect_equal(expected_gr)        
} )
