library(testthat)

library(tidyverse)
library(plyranges)

source(test_path('../../StemCNV-check/scripts/R/R_io_functions.R'))

#  Functions to test:
# - [ ] get_chromosome_set
# - [ ] get_sample_info
# - [ ] read_raw     #! replace with read_SNP_vcf at some point
# - [ ] read_PennCNV #! replace with read_CNV_vcf at some point
# - [ ] read_CBS     #! replace with read_CNV_vcf at some point
# - [ ] load_preprocessed_cnvs #! consider replacing this with read_CNV_vcf (+ annotation support)
# - [ ] get_static_path #! hard to test; should only do things if insside the Rmd environment
# - [x] load_gtf_data
# - [x] load_genomeInfo
# - [ ] ensure_list_cols
# - [ ] get_SNP_clustering_IDs #! this should go to a different script

# add test with completely empty callset
test_that("load_gtf_data", {
  # Construct minimal config
  config <- list(
    'static_data' = list(
      'genome_gtf_file' = test_path('../data/hg_minimal.gtf')
    ),
    'settings' = list(
      'CNV_processing' = list(
        'gene_overlap' = list(
          'exclude_gene_type_regex' = c(),
          'include_only_these_gene_types' = c('lncRNA', 'miRNA', 'protein_coding')
        )
      )
    )
  )
  expected_tb <- tibble(
    seqnames = 'chr1',
    start = c(4000, 11873, 20000, 28050000, 50900000),
    end = c(5000, 14409, 50000, 28070000, 50910000),
    strand = '+',
    source = factor(NA, levels = character()),
    type = factor('gene', levels = c('gene', 'transcript', 'exon')),
    gene_id = c("dummy_1", 'ENSG00000223972', 'dummy_2' , 'dummy_3', 'dummy_4'),
    gene_type = c("unprocessed_pseudogene", 'lncRNA', 'something_else', 'protein_coding', 'protein_coding'),
    gene_name = c("dummyA", 'DDX11L1', 'dummyB', 'dummyC', 'dummyD')
  )

  expect_equal(load_gtf_data(config), as_granges(expected_tb %>% filter(gene_type %!in% c("unprocessed_pseudogene", 'something_else'))))
  config$settings$CNV_processing$gene_overlap$include_only_these_gene_types <- c('lncRNA', 'miRNA', 'protein_coding', 'something_else')
  config$settings$CNV_processing$gene_overlap$exclude_gene_type_regex <- c('lncRNA')
  expect_equal(load_gtf_data(config), as_granges(expected_tb %>% filter(gene_type %!in% c('lncRNA', 'unprocessed_pseudogene'))))
  config$settings$CNV_processing$gene_overlap$include_only_these_gene_types <- c()
  config$settings$CNV_processing$gene_overlap$exclude_gene_type_regex <- c('lncRNA', 'something_else', 'unprocessed_pseudogene')
  expect_equal(load_gtf_data(config), as_granges(expected_tb %>% filter(gene_type == 'protein_coding')))
  config$settings$CNV_processing$gene_overlap$include_only_these_gene_types <- c()
  config$settings$CNV_processing$gene_overlap$exclude_gene_type_regex <- c()
  expect_equal(load_gtf_data(config), as_granges(expected_tb))
} )

test_that('load_genomeInfo', {
  config <- list(
    'static_data' = list(
      'genomeInfo_file' = test_path('../data/gr_info_minimal.tsv')
    )
  )
  expected_gr <- tibble(
    seqnames = 'chr1',
    size = 249250621,
    start = c(0, 2300000, 5400000, 28000000, 30200000, 32400000, 46800000, 120600000, 121500000, 125000000, 128900000, 142600000, 147000000, 150300000),
    end = c(2300000, 5400000, 7200000, 30200000, 32400000, 34600000, 50700000, 121500000, 125000000, 128900000, 142600000, 147000000, 150300000, 155000000),
    band_name = c("p36.33", "p36.32", "p36.31", "p35.3", "p35.2", "p35.1", "p33", "p11.2", "p11.1", "q11", "q12", "q21.1", "q21.2", "q21.3"),
    band_staining = c("gneg", "gpos25", "gneg", "gpos25", "gneg", "gpos25", "gpos75", "gneg", "acen", "acen", "gvar", "gneg", "gpos50", "gneg"),
    centromer = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE),
    section_name = c("1p36.33", "1p36.32", "1p36.31", "1p35.3", "1p35.2", "1p35.1", "1p33", "1p11.2", "1p11.1", "1q11", "1q12", "1q21.1", "1q21.2", "1q21.3")
  ) %>% as_granges()

  expect_equal(load_genomeInfo(config), expected_gr)
})