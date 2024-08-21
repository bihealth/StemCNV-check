library(testthat)

library(tidyverse)
library(plyranges)
library(vcfR)

source(test_path('../../stemcnv_check/scripts/R/vcf_io_functions.R')) # for vcfR_to_tibble
source(test_path('../../stemcnv_check/scripts/R/preprocess_CNV_functions.R'))

# Example calls:
# - <1000bp,    will be filtered
# - <5 probes, will be filtered
# - <10 probes/kb, will be filtered
# - matching pair within 500bp, will be merged
# - pair within 500bp w/o CNV_type match, will not be merged
# - pair that will be filtered, after merging due to density
# - pair of calls <1000bp, but will be merged and kept then
# - 3 calls that will be kept & merged

snp.vcf <- parse_snp_vcf(test_path('../data/minimal_probes.vcf')) 

# minimal / expected cols:
# seqnames, start, end, width, ID, CNV_caller, CNV_type, CN, sample_id

cnv_tb_raw <- tibble(
  sample_id = 'test_sample',
  seqnames  = rep('chr1', 14),
  # Note: granges are 1-based & inclusive, so end is -1 to have matching width
  start = c(100, 1000, 1.0e8, 3000, 4400, 6000, 7400, 1e6, 2e6+400, 9000,  9600, 12e4, 13e4, 14e4) %>% as.integer(),
  end   = c(200, 1600, 1.1e8, 4000, 5400, 7000, 8400, 2e6, 3e6+400, 9400, 10000, 13e4, 14e4, 15e4)-1 %>% as.integer(),
  width = c(100, 600,    1e7, 1000, 1000, 1000, 1000, 1e6, 1e6,      400,   400,  1e4,  1e4,  1e4),
  CNV_caller = 'Test',
  CNV_type = c(rep('DUP', 6), rep('DEL', 8)),
  CN = c(rep(3, 6),rep(1, 5), 1, 0, 1),
  ID = paste(CNV_caller, CNV_type, seqnames, start, end, sep='_'),
)

cnv_tb_snps <- cnv_tb_raw
cnv_tb_snps$n_probes <-      c(5, 3, 20, 6, 5, 5, 5, 10, 10, 5, 5, 5, 5, 5)
cnv_tb_snps$n_uniq_probes <- c(5, 3, 20, 5, 5, 5, 5, 10, 10, 5, 5, 5, 5, 5)
cnv_tb_snps$probe_density_Mb <- cnv_tb_snps$n_uniq_probes / cnv_tb_snps$width * 1e6


merged_tb <- tibble(
  seqnames  = rep('chr1', 9),
  start = c(100, 1000, 3000, 6000, 1.0e8, 7400,  9000, 12e4, 1e6) %>% as.integer(),
  end   = c(200, 1600, 5400, 7000, 1.1e8, 8400, 10000, 15e4, 3e6+400)-1 %>% as.integer(),
  width = c(100, 600, 2400, 1000, 1e7, 1000, 1000, 3e4, 2e6+400),
  sample_id = 'test_sample',
  CNV_caller = 'Test',
  CNV_type = c(rep('DUP', 5), rep('DEL', 4)),
  n_initial_calls = c(1, 1, 2, 1, 1, 1, 2, 3, 2),
  initial_call_details = c(NA, NA, '3000_3999_CN3;4400_5399_CN3', NA, NA, 
                           NA, '9000_9399_CN1;9600_9999_CN1', '120000_129999_CN1;130000_139999_CN0;140000_149999_CN1',
                            '1000000_1999999_CN1;2000400_3000399_CN1'),
  CN = c(rep(3, 5), rep(1, 4)),  
  ID = paste(CNV_caller, CNV_type, seqnames, start, end, sep='_'),
  n_probes =      c(5, 3, 11, 5, 20, 5, 10, 15, 20),
  n_uniq_probes = c(5, 3, 10, 5, 20, 5, 10, 15, 20),
  probe_density_Mb = n_uniq_probes / width * 1e6,
  
)

merged_filtered_tb <- merged_tb 
merged_filtered_tb$FILTER <- c('Size', 'Size;n_probes', 'PASS', 'PASS',
                               'Density', 'PASS', 'PASS', 'PASS', 'Density')

# Functions to test
# - add_snp_probe_counts (only gr, + snp_vcf)
# - merge_calls (accepts tb or gr), runs add_snp_probe_counts
# - add_call_prefilters (only gr)
# - apply_preprocessing (tb or gr, + snp_vcf), runs merge_calls, add_call_prefilters
# - get_median_LRR (only gr, + snp_vcf)

# minimal_probes.vcf test file:
# - expected number of probes in calls
# - 1 call with FILTER set
# - 1 pos with dupl. position
# - 2 probes outside calls
# - Future: add actual LRR values to test vcf writing?

# Note: tests won't catch issues with mismatched CHROM names
test_that("add snp probe counts", {
    as_granges(cnv_tb_raw) %>%
        add_snp_probe_counts(snp.vcf) %>%
        expect_equal(as_granges(cnv_tb_snps))    
})


test_that("merge & add snp counts calls", {
    merge.distance <- 500
    merge_calls(cnv_tb_raw, merge.distance, snp.vcf) %>%
        expect_equal(as_granges(merged_tb))
})

test_that("merge then filter calls", {
    tool_config <- list(
        filter.minlength = 1000,
        filter.minprobes = 5,
        filter.mindensity.Mb = 10
    )
    merge.distance <- 500
    merge_calls(cnv_tb_raw, merge.distance, snp.vcf) %>%
        add_call_prefilters(tool_config) %>%
        expect_equal(as_granges(merged_filtered_tb))  
} )

# Deprecated, (pre)filtering now always happens after merging
#
# test_that("filter raw calls", { } )
# test_that("filter then merge calls", { } )

test_that("test empty calls", {
    tool_config <- list(
        filter.minlength = 1000,
        filter.minprobes = 5,
        filter.mindensity.Mb = 10
    )
    merge.distance <- 500
    
    empty_raw <- cnv_tb_raw %>%
        filter(sample_id == 'non_existent_sample')
    empty_snps <- cnv_tb_snps %>%
        filter(sample_id == 'non_existent_sample')
    empty_merged <- merged_tb %>%
        filter(sample_id == 'non_existent_sample') %>%
        as_granges()
    
    as_granges(empty_raw) %>%
        add_snp_probe_counts(snp.vcf) %>%
        expect_equal(as_granges(empty_snps))
    
    merge_calls(empty_raw, merge.distance, snp.vcf) %>%
        expect_equal(empty_merged)
    
    empty_merged$FILTER <- character()
    merge_calls(empty_raw, merge.distance, snp.vcf) %>%
        add_call_prefilters(tool_config) %>%
        expect_equal(as_granges(empty_merged))  
    
} )

# Test fix_CHROM_format (from helper functions) with these example tbs
test_that('fix_CHROM_format', {
    levels <- c(1:5,7,20,22,'X','Y','MT') %>% as.character()
    tb <- cnv_tb_raw %>%
        mutate(
            seqnames = c(1,2,4,'X','X','Y',20,2,7,22,5,3,'MT','X') %>% 
                as.character() #%>%
                #factor(levels = levels)
        )
    expected_tb <- tb %>%
        mutate(seqnames = factor(seqnames, levels = levels))

    expect_equal(
        fix_CHROM_format(as_granges(tb), 'NCBI'), 
        expected_tb %>% as_granges()
    )
    expect_equal(
        fix_CHROM_format(as_granges(tb), 'Ensembl'),
        expected_tb %>% as_granges()
    )
    
    expected_tb <- tb %>%
        mutate(
            seqnames = paste0('chr', str_remove(seqnames, 'T')) %>% 
                factor(levels = paste0('chr', str_remove(levels, 'T')))
        )    
    expect_equal(
        fix_CHROM_format(as_granges(tb), 'UCSC'), 
        expected_tb %>% as_granges()
    )
    
    expect_error(fix_CHROM_format(as_granges(tb), 'unknown_style'))
})


# apply_preprocessing
test_that("apply_preprocessing", {
    
    tool_config <- list(
        merge.distance = 500,
        filter.minlength = 1000,
        filter.minprobes = 5,
        filter.mindensity.Mb = 10
    )
    
    expected <- merged_filtered_tb %>%
        mutate(
            seqnames = str_remove(seqnames, 'chr') %>% factor(),
            ID = paste(CNV_caller, CNV_type, seqnames, start, end, sep='_')
        ) %>%
        arrange(seqnames, start) %>%
        as_granges()

    # Mismatched CHROM style
    expect_error(
        as_granges(cnv_tb_raw) %>%
            fix_CHROM_format('NCBI') %>%
            apply_preprocessing(
                snp.vcf, 
                tool_config
            )
    )
   
    
    snp.vcf <- fix_CHROM_format(snp.vcf, 'NCBI')
    as_granges(cnv_tb_raw) %>%
            fix_CHROM_format('NCBI') %>%
            apply_preprocessing(
                snp.vcf, 
                tool_config
            ) %>% 
        expect_equal(expected)
    
})

test_that('get_median_LRR', {
    
    expected_LRR <- merged_tb %>%
        mutate(LRR = c(1, 1.3, 0.89, 2, 1.385, 0, -0.88, -1.31, -0.895))
    
    merged_tb %>%
        as_granges() %>%
        get_median_LRR(snp.vcf) %>%
        expect_equal(as_granges(expected_LRR))
    
})

# test BAF cluster determination
# Need to add cases where this should fail
#BAF clusters:
#4, 2, 4, 3, 4, 2, 2, 2, 2  