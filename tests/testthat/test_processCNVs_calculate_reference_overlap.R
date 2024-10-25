library(testthat)

library(tidyverse)
library(plyranges)

source(test_path('../../stemcnv_check/scripts/R/helper_functions.R'))
source(test_path('../../stemcnv_check/scripts/R/processCNVs_calculate_reference_overlap.R'))

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
  CNV_caller = c('StemCNV-check', 'toolA', 'StemCNV-check', 'StemCNV-check', 'toolA', 'toolB'),
  n_probes = c(15, 5, 15, 25, 5, 5),
  CN = c(3, 3, 0, 1, 1, 1), 
)

ref_cnvs <- tibble(   
  seqnames = 'chr1',
  start = c(4000, 10000, 12000, 18e4, 25e4, 45e4) %>% as.integer(),
  end   = c(5500, 10850, 14000, 20e4, 30e4, 50e4) %>% as.integer(),
  sample_id = 'test_sample',
  CNV_type = c('DUP', 'DUP', 'DUP', 'DEL', 'DEL', 'DEL'),
  ID = paste('combined', CNV_type, seqnames, start, end, sep='_'),
  CNV_caller = c('StemCNV-check', 'faketool', 'StemCNV-check', 'StemCNV-check', 'toolA', 'toolB'),
  n_probes = c(15, 20, 15, 25, 5, 5),
  CN = c(3, 4, 3, 1, 1, 1),
) %>%
  as_granges()

test_that("annotate_reference_overlap", {
  min.reciprocal.coverage.with.ref <- 80
  
  expected_gr <- sample_cnvs %>%
    mutate(
        reference_overlap = c(T, T, F, F, F, F),
        reference_coverage = c(100, 85.01, NA_real_, NA_real_, NA_real_, NA_real_),
        reference_caller = c('StemCNV-check', 'faketool', NA_character_, NA_character_,NA_character_,NA_character_)
    ) %>%
    as_granges()
  
  annotate_reference_overlap(
      as_granges(sample_cnvs),
      ref_cnvs,
      min.reciprocal.coverage.with.ref
  ) %>% 
    expect_equal(expected_gr) 
} )

test_that("annotate_reference_overlap, empty ref", {
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

test_that("annotate_reference_overlap, empty input", {
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
