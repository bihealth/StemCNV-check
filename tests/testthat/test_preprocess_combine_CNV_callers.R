library(testthat)

library(tidyverse)
library(plyranges)


source(test_path('../../StemCNV-check/scripts/R/R_io_functions.R'))
source(test_path('../../StemCNV-check/scripts/R/preprocess_combine_CNV_callers.R'))

# Example calls (toolA and toolB):

# 0+1 - unqiue in one tool
# 1+0 - unqiue in one tool
# 1+1 - same call in both tools
# 1+1 - same region, but diff CNV_type
# 1+1 - both tools, but median region cov <60% (100%: 2k, <20% : <400)
# 1+1 - same call, but different copynumber

# 2+1 - 3 calls, all merged
# 2+1 - 3 calls, region with 2 calls; median region cov <60%  (100%: 9k, <20% : <1.8k)
# 2+1 - 3 calls, region with 1 call ; single call <50% region cov (4k---4k + 3k)
# 2+2 - 4 calls, both regions >60% median cov, but no single call >50% cov

toolA <- tibble(
  seqnames = rep('chr1', 13),
  start = c(1000, 4000, 6000,  8000, 12000, 15e4, 18e4, 22.2e4, 27e4, 31e4, 36e4, 50e4, 55e4) %>% as.integer(),
  end   = c(2000, 5500, 7500, 10000, 14000, 17e4, 20e4, 23e4, 27.8e4, 35e4, 40e4, 54e4, 58e4) %>% as.integer(),
  CNV_type = c(rep('gain', 3), rep('loss', 10)),
  sample_id = 'test_sample',
  CNV_caller = 'toolA',
  n_premerged_calls = 2,
  n_snp_probes = 10,
  copynumber = c(rep('3', 3), rep('1', 10)) %>% as.character(),
  caller_confidence = 1 ,
  ID = paste(CNV_caller, CNV_type, seqnames, start, end, sep='_'),
)

toolB <- tibble(
  seqnames = rep('chr1', 10),
  start = c(2000, 4000, 6000,  9700, 12000, 15e4, 21e4, 33e4, 52e4, 57e4) %>% as.integer(),
  end   = c(3000, 5500, 7500, 10000, 14000, 20e4, 30e4, 36e4, 56e4, 60e4) %>% as.integer(),
  CNV_type = c(rep('gain', 2), rep('loss', 8)),
  sample_id = 'test_sample',
  CNV_caller = 'toolB',
  n_premerged_calls = 1,
  n_snp_probes = 5,
  copynumber = c(rep('3', 2), '1', '1', '0', rep('1', 5)) %>% as.character(),
  caller_confidence = 1 ,
  ID = paste(CNV_caller, CNV_type, seqnames, start, end, sep='_'),
)


combined_tools <- bind_ranges(
  #combined calls
  tibble(   
    seqnames = 'chr1',
    start = c(4000, 12000, 15e4) %>% as.integer(),
    end   = c(5500, 14000, 20e4) %>% as.integer(),
    sample_id = 'test_sample',
    CNV_type = c('gain', 'loss', 'loss'),
    ID = paste('combined', CNV_type, seqnames, start, end, sep='_'),
    CNV_caller = list(c('toolA','toolB'), c('toolA','toolB'), c('toolA', 'toolA','toolB')),
    n_premerged_calls = list(c(2,1), c(2,1), c(2,2,1)),
    n_snp_probes = list(c(10,5), c(10,5), c(10,10,5)),
    copynumber = list(c('3','3'), c('1','0'), c('1','1','1')),
    caller_confidence = list(c(1,1), c(1,1), c(1,1,1)),
    overlap_merged_call = NA_real_,
    caller_merging_coverage = c('toolA-100,toolB-100', 'toolA-100,toolB-100', 'toolA-80,toolB-100'),
    caller_merging_state = 'combined'
  ) %>% as_granges(),
  # single calls
  # sorted by type, then caller, then start
  bind_rows(
    toolA[c(1,3,4,8:13),],
    toolB[c(1,3,4,7:10),]
  ) %>%
    arrange(CNV_type, CNV_caller, start) %>%
    mutate(overlap_merged_call = NA_real_,
           caller_merging_coverage = NA_character_,
           caller_merging_state = 'no-overlap') %>%
    ensure_list_cols(),
  # pre.ov calls
  bind_rows(
    toolA[c(2,5,6,7),],
    toolB[c(2,5,6),]
  ) %>%
    arrange(CNV_type, CNV_caller, start) %>%
    mutate(
      # gain-A-4500, gain-B-4500, loss-A-12000, loss-A-15000, loss-A-18000, loss-B-12000, loss-B-15000
      overlap_merged_call = 100*c(1, 1, 1, (20001)/(50001), (20001)/(50001), 1, 1) %>% as.double(),
      caller_merging_coverage = NA_character_,
      caller_merging_state = 'pre-overlap') %>%
    ensure_list_cols()
)

test_that("combined 2 CNV callers", {
  min.greatest.region.overlap <- 50
  min.median.tool.coverage    <- 60
  
  list(
    'toolA' = as_granges(toolA),
    'toolB' = as_granges(toolB)
  ) %>%
    bind_ranges() %>%
    combine_CNV_callers(min.greatest.region.overlap, min.median.tool.coverage) %>%
    expect_equal(combined_tools)        
} )
