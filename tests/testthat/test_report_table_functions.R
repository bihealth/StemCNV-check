library(tidyverse)
library(plyranges)
library(testthat)

source(test_path("../../StemCNV-check/scripts/R/R_io_functions.R"))
source(test_path("../../StemCNV-check/scripts/R/R_table_functions.R"))

# Test `format_hotspots_to_badge` function
# format_hotspots_to_badge <- function(hotspot_vec, CNVtype_vec, gene_details, listname = 'high_impact')
test_that("format_hotspots_to_badge", {
  local_edition(3)
  input_hotspot_vector <- c("", "12p13.3", "12p13.3", "TP53", "18q21,JAK2", "18q21,JAK2")
  input_CNVtype_vector <- c("gain", "gain", "LOH", "gain", "loss", "gain")
  # 1 - empty
  # 2 - gband hit, matching CNV (gain)
  # 3 - gband hit, not matching CNV (Note: this case should never occur; warning needed?)
  # 4 - gene hit
  # 5 - gene hit & gband hit matching CNV (loss)
  # 6 - gene hit & gband hit not matching CNV (Note: theoretically possible)

  gene_details <- read_tsv(test_path('../../StemCNV-check/supplemental-files/HighImpact-stemcell-hotspots.tsv'))
  
  # Compare the output with the expected output
  expect_snapshot(cat(paste(
    format_hotspots_to_badge(input_hotspot_vector, input_CNVtype_vector, gene_details, 'high_impact'),
    collapse = "\n"
  )))
})