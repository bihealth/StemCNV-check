library(tidyverse)
library(plyranges)
library(testthat)

source(test_path("../../stemcnv_check/scripts/R/helper_functions.R"))
source(test_path("../../stemcnv_check/scripts/R/R_table_functions.R"))

# Functions to test:
# - [ ] vector_to_js
# - [ ] format_column_names
# - [ ] simple_table_output
# - [ ] summary_table
# - [x] format_hotspots_to_badge
# - [ ] CNV_table_output
# - [ ] gene_table_output

# Test `format_hotspots_to_badge` function
# format_hotspots_to_badge <- function(hotspot_vec, CNVtype_vec, gene_details, listname = 'high_impact')
test_that("format_hotspots_to_badge", {
  testthat::local_edition(3)
  input_hotspot_vector <- c("", "12p13.3", "12p13.3", "TP53", "18q21,JAK2", "18q21,JAK2")
  input_CNVtype_vector <- c("gain", "gain", "LOH", "gain", "loss", "gain")
  # 1 - empty
  # 2 - gband hit, matching CNV (gain)
  # 3 - gband hit, not matching CNV (Note: this case should never occur; warning needed?)
  # 4 - gene hit
  # 5 - gene hit & gband hit matching CNV (loss)
  # 6 - gene hit & gband hit not matching CNV (Note: theoretically possible)

  gene_details <- read_tsv(test_path('../../stemcnv_check/supplemental-files/HighImpact-stemcell-hotspots.tsv'), 
                           show_col_types = FALSE)
  
  expected <- c(
    '-', 
    '<span class="badge badge-HI" title="high_impact list name: StemCell-Hotspots&#013;Annotation source:&#013;StemCNV-check curation; ISCCR guidelines 2023&#013; doi.org/10.1038/ncomms5825 &#013; https://doi.org/10.1101/2021.05.22.445238">12p13.3</span>', 
    '12p13.3', 
    '<span class="badge badge-HI" title="high_impact list name: StemCell-Hotspots&#013;Annotation source:&#013;StemCNV-check curation; doi.org/10.1038/s41588-022-01147-3">TP53</span>', 
    '<span class="badge badge-HI" title="high_impact list name: StemCell-Hotspots&#013;Annotation source:&#013;StemCNV-check curation; ISCCR guidelines 2023">18q21</span><span class="badge badge-HI" title="high_impact list name: StemCell-Hotspots&#013;Annotation source:&#013;StemCNV-check curation; doi.org/10.1038/s41588-022-01147-3">JAK2</span>', 
    '18q21<span class="badge badge-HI" title="high_impact list name: StemCell-Hotspots&#013;Annotation source:&#013;StemCNV-check curation; doi.org/10.1038/s41588-022-01147-3">JAK2</span>'
  ) 
  
  expect_equal(
    format_hotspots_to_badge(input_hotspot_vector, input_CNVtype_vector, gene_details, 'high_impact'),
    expected
  )
})