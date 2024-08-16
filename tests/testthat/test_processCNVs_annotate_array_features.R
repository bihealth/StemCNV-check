library(testthat)

library(tidyverse)
library(plyranges)

source(test_path('../../stemcnv_check/scripts/R/helper_functions.R'))
source(test_path('../../stemcnv_check/scripts/R/processCNVs_annotate_array_features.R'))

# Functions:
# - get_accurate_snp_probe_count !Note: this will be phased out with vcf overhauls
# - annotate_gaps
# - annotate_high_density

sample_cnvs <- tibble(
  seqnames = 'chr1',
  start = c(4000, 10000, 28000000, 28060000, 40000, 5000000, 31000000),
  end   = c(5500, 14000, 28055000, 28065000, 50000, 7000000, 32000000),
  sample_id = 'test_sample',
  CNV_type = c('gain', 'gain', 'gain', 'loss', 'loss', 'loss', 'loss'),
  ID = paste('combined', CNV_type, seqnames, start, end, sep='_'),
  CNV_caller = list(c('toolA','toolB'), 'toolA', c('toolA','toolB'), c('toolA','toolA','toolB'), 'toolA', 'toolB', 'toolB'),
  # n_premerged_calls = list(c(2,1), 1, c(2,1), c(2,2,1),1,1,1),
  n_probes = c(15, 100, 150, 250, 100, 50, 50),
  n_uniq_probes = c(15, 100, 150, 100, 100, 50, 50),
  CN = c(3, 3, 4, 1, 1, 1, 0),
  # caller_confidence = list(c(1,1), 1, c(1,1), c(1,1,1), 1,1,1),
  overlap_merged_call = NA_real_,
  caller_merging_coverage = c('toolA-100,toolB-100', NA, 'toolA-100,toolB-100', 'toolA-80,toolB-100', NA, NA, NA),
  caller_merging_state = c('combined', NA, rep('combined', 2), 'no-overlap', 'no-overlap', 'no-overlap'),
  reference_overlap = c(T, T, F, F, F, F, T),
  reference_coverage = c(100, 85.01, NA_real_, NA_real_, NA_real_, NA_real_, 60) %>% as.list(),
  reference_caller = list(c('toolA','toolB'), 'faketool', NA_character_, NA_character_,NA_character_,NA_character_, 'toolA'),
  test_hits = c(NA, 'DDX11L1', 'dummyC', NA, '1p36,chr1:40000-50000', '1p36', '1p35.2')
) %>% as_granges()

# Gaps:
# - no gaps in CNV
# - single gap in CNV         -> above th, true
# - multiple gaps in one CNV  -> above th, true
# - CNV has single gap,    but below min.perc.gap_area
# - CNV has multiple gaps, but below min.perc.gap_area
# - CNV has gaps, but ! (gap_slope * percent_gap_coverage + gap_intercept) <= log2(n_uniq_probes)

# - gap overlaps CNV border -> error
# - gap fully conatins CNV  -> error


test_that("Annotate CNVs with gaps", {
  gapfile <- test_path('../data/gaps_minimal.bed')
  gap_area.uniq_probes.rel <- list(-12, 12.5)
  min.perc.gap_area <- 0.33
  
  expexted_gr <- sample_cnvs %>%
    mutate(
      percent_gap_coverage = c(0, 2000/4001, 25000/55001, 1000/5001, 2000/10001, 1e6/(2e6+1), 0),
      probe_coverage_gap = c(FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE)
    )
  
  annotate_gaps(sample_cnvs, gapfile, min.perc.gap_area, gap_area.uniq_probes.rel) %>%
    expect_equal(expexted_gr)        
  
  # Call partially overlapping a gap (one endpoint in gap)
  expect_error(
    annotate_gaps(
      sample_cnvs %>% bind_ranges(
            tibble(seqnames = 'chr1', start = 12000, end = 15000, ID = 'error_test') %>%
              as_granges()
      ),
      gapfile, min.perc.gap_area, gap_area.uniq_probes.rel
    ),
    regexp = 'CNV call endpoint\\(s\\) overlap with gap areas from ".+/data/gaps_minimal.bed": error_test'
  )
  
} )

test_that("Annotate CNVs with probe density flags", {
  density_file <- test_path('../data/density_minimal.bed')
  density.quantile.cutoff <- 0.99
  
  expexted_gr <- sample_cnvs %>%
    mutate(
      high_probe_density = c(NA, NA, TRUE, TRUE, NA, FALSE, FALSE)
    )
  
  annotate_high_density(sample_cnvs, density_file, density.quantile.cutoff) %>%
    expect_equal(expexted_gr)        

} )