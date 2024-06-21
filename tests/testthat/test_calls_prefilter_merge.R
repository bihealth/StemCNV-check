library(testthat)

library(tidyverse)
library(plyranges)

source('StemCNV-check/scripts/R/calls_prefilter_merge.R')

test_that("prefilter a tibble", {
  merge.distance <- 500
  # setup tibble with mock calls; ensure some calls are within merge.distance, some not
  df <- tibble(
      CNV_type = c('gain', 'gain', 'gain', 'loss'),
      Chr = c('chr1', 'chr1', 'chr1', 'chr1'),
      start = c(100, 200, 1000, 2000),
      end = c(200, 300, 1100, 2100),
      CNV_caller = 'Test',
      n_snp_probes = c(5, 5, 5, 5),
      copynumber = c(3, 3, 3, 1),
      caller_confidence = c(0.9, 0.9, 0.9, 0.9),
  )

  merge_calls(df, merge.distance) %>%
    as_tibble() %>%
    expect_equal(
      tibble(
        CNV_type = c('gain', 'loss', 'gain', 'loss'),
        Chr = c('chr1', 'chr1', 'chr1', 'chr1'),
        start = c(100, 1000, 2000),
        end = c(300, 1100, 2100),
        CNV_caller = 'Test',
        ID = c('Test_gain_chr1_100_300', 'Test_gain_chr1_1000_1100', 'Test_loss_chr1_2000_2100'),
        n_premerged_calls = c(2, 1, 1),
        n_snp_probes = c(10, 5, 5),
        copynumber = c(3, 3, 1),
        caller_confidence = c(0.9, 0.9, 0.9)
      )
    )


})