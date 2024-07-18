library(testthat)

library(tidyverse)
library(plyranges)

source(test_path('../../stemcnv_check/scripts/R/processCNVs_calls_prefilter_merge.R'))

# Example calls:
# - <1000bp,    will be filtered
# - <5 probes, will be filtered
# - <10 probes/kb, will be filtered
# - matching pair within 500bp, will be merged
# - pair within 500bp w/o CNV_type match, will not be merged
# - pair that will be filtered, after merging due to density
# - pair of calls <1000bp, but will be merged and kept then
# - 3 calls that will be kept & merged

raw_tb <- tibble(
  sample_id = 'test_sample',
  Chr   = rep('chr1', 14),
  start = c(100, 1000, 1.0e8, 3000, 4400, 6000, 7400, 1e6, 2e6+500, 9000,  9600, 12e4, 13e4, 14e4),
  end   = c(200, 1600, 1.1e8, 4000, 5400, 7000, 8400, 2e6, 3e6+500, 9400, 10000, 13e4, 14e4, 15e4),
  length = c(100, 600,   1e7, 1000, 1000, 1000, 1000, 1e6, 1e6,      400,   400,  1e4,  1e4,  1e4),
  n_snp_probes = c(5, 3,  20,    5,    5,    5,    5,   10, 10,        5,     5,    5,    5,    5),
  snp.density = n_snp_probes / length * 1e6,
  CNV_type = c(rep('gain', 6), rep('loss', 8)),
  copynumber = c(rep(3, 6), rep(1, 8)),
  CNV_caller = 'Test',
  ID = paste(CNV_caller, CNV_type, Chr, start, end, sep='_'),
  caller_confidence = 1
)

filtered_tb <- tibble(
  sample_id = 'test_sample',
  Chr   = rep('chr1', 9),
  start = c(3000, 4400, 6000, 7400, 1e6, 2e6+500, 12e4, 13e4, 14e4),
  end   = c(4000, 5400, 7000, 8400, 2e6, 3e6+500, 13e4, 14e4, 15e4),
  length = c(1000,1000, 1000, 1000, 1e6, 1e6,      1e4,  1e4,  1e4),
  n_snp_probes = c(rep(5, 4), 10, 10, 5,5,5),
  snp.density = n_snp_probes / length * 1e6,
  CNV_type = c(rep('gain', 3), rep('loss', 6)),
  copynumber = c(rep(3, 3), rep(1, 6)),
  CNV_caller = 'Test',
  ID = paste(CNV_caller, CNV_type, Chr, start, end, sep='_'),
  caller_confidence = 1
)

merged_tb <- tibble(
  Chr   = rep('chr1', 9),
  start = c(100, 1000, 3000, 6000, 1.0e8, 7400,  9000, 12e4, 1e6) %>% as.integer(),
  end   = c(200, 1600, 5400, 7000, 1.1e8, 8400, 10000, 15e4, 3e6+500) %>% as.integer(),
  CNV_type = c(rep('gain', 5), rep('loss', 4)),
  sample_id = 'test_sample',
  CNV_caller = 'Test',
  #length = c(100, 600, 2400, 1000,   1e7, 1000, 2e6+500) - 1,
  n_premerged_calls = c(1, 1, 2, 1, 1, 1, 2, 3, 2),
  n_snp_probes = c(5, 3, 10, 5, 20, 5, 10, 15, 20),
  copynumber = c(rep(3, 5), rep(1, 4)) %>% as.character(),
  caller_confidence = 1 ,
  ID = paste(CNV_caller, CNV_type, Chr, start, end, sep='_'),
)

filtered_merged_tb <- tibble(
  Chr   = rep('chr1', 5),
  start = c(3000, 6000, 7400, 12e4, 1e6) %>% as.integer(),
  end   = c(5400, 7000, 8400, 15e4, 3e6+500) %>% as.integer(),
  CNV_type = c('gain', 'gain', 'loss', 'loss', 'loss'),
  sample_id = 'test_sample',
  CNV_caller = 'Test',
  #length = c(2000,1000, 1000, 3000, 2e6+500),
  n_premerged_calls = c(2, 1, 1, 3, 2),
  n_snp_probes = c(10, 5, 5, 15, 20),  
  copynumber = c(rep(3, 2), rep(1, 3))  %>% as.character(),
  caller_confidence = 1,
  ID = paste(CNV_caller, CNV_type, Chr, start, end, sep='_'),
)

merged_filtered_tb <- tibble(
  Chr   = rep('chr1', 5),
  start = c(3000, 6000, 7400,  9000, 12e4) %>% as.integer(),
  end   = c(5400, 7000, 8400, 10000, 15e4) %>% as.integer(),
  CNV_type = c(rep('gain', 2), rep('loss', 3)),
  sample_id = 'test_sample',
  CNV_caller = 'Test',
  #length = c(2400, 1000, 1000, 1000, 3e4) - 1,
  n_premerged_calls = c(2, 1, 1, 2, 3),
  n_snp_probes = c(10, 5, 5, 10, 15),
  copynumber = c(rep(3, 2), rep(1, 3)) %>% as.character(),
  caller_confidence = 1 ,
  ID = paste(CNV_caller, CNV_type, Chr, start, end, sep='_'),
)

test_that("filter raw calls", {
  min.length <- 1000
  min.snps <- 5
  min.dens <- 10
  prefilter_calls(raw_tb, min.snps, min.length, min.dens) %>%
    expect_equal(filtered_tb)        
} )

test_that("merge raw calls", {
  merge.distance <- 500
  merge_calls(raw_tb, merge.distance) %>%
    expect_equal(as_granges(merged_tb, seqnames=Chr))
})

test_that("filter then merge calls", {
  min.length <- 1000
  min.snps <- 5
  min.dens <- 10
  merge.distance <- 500
  prefilter_calls(raw_tb, min.snps, min.length, min.dens) %>%
    merge_calls(merge.distance) %>%
    expect_equal(as_granges(filtered_merged_tb, seqnames=Chr))  
} )

test_that("merge then filter calls", {
  min.length <- 1000
  min.snps <- 5
  min.dens <- 10
  merge.distance <- 500
  merge_calls(raw_tb, merge.distance) %>%
    prefilter_calls(min.snps, min.length, min.dens) %>%
    expect_equal(as_granges(merged_filtered_tb, seqnames=Chr))  
} )

# add tests for empty callsets
empty_tb <- raw_tb %>%
  filter(sample_id == 'non_existent_sample')
empty_gr <- merged_tb %>%
  dplyr::rename(seqnames = Chr) %>%
  filter(sample_id == 'non_existent_sample') %>%
  as_granges()

test_that("filter empty calls", {
  min.length <- 1000
  min.snps <- 5
  min.dens <- 10
  
  prefilter_calls(empty_tb, min.snps, min.length, min.dens) %>%
    expect_equal(empty_tb) 
  
  prefilter_calls(empty_gr, min.snps, min.length, min.dens) %>%
    expect_equal(empty_gr)
} )

test_that("merge empty calls", {
  merge.distance <- 500

  merge_calls(empty_tb, merge.distance) %>%
    expect_equal(empty_gr)        
  
  merge_calls(empty_gr, merge.distance) %>%
    expect_equal(empty_gr)
} )