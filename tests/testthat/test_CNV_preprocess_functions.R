library(testthat)

library(tidyverse)
library(plyranges)
library(vcfR)

source(test_path('../../stemcnv_check/scripts/R/vcf_io_functions.R')) # for vcfR_to_tibble
source(test_path('../../stemcnv_check/scripts/R/CNV_preprocess_functions.R'))

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

# Example calls:
# - <1000bp,    will be filtered
# - <5 probes, will be filtered
# - <10 probes/kb, will be filtered
# - matching pair within 500bp, will be merged
# - pair within 500bp w/o CNV_type match, will not be merged
# - pair that will be filtered, after merging due to density
# - pair of calls <1000bp, but will be merged and kept then
# - 3 calls that will be kept & merged

snp_vcf_gr <- parse_snp_vcf(test_path('../data/minimal_probes.vcf')) 

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
merged_filtered_tb$FILTER <- c('min_size', 'min_size;min_probes', NA, NA,
                               'min_density', NA, NA, NA, 'min_density')

test_that("add_snp_probe_counts", {
    as_granges(cnv_tb_raw) %>%
        add_snp_probe_counts(snp_vcf_gr) %>%
        expect_equal(as_granges(cnv_tb_snps))    
})


test_that("merge_calls", {
    # test absolute distance only
    merge_config <- list(
        maximum.gap.allowed = 1000,
        call.extension.percent = 0,
        merge.gap.absolute = 500,
        merge.gap.snps = 0
    )
    merge_calls(cnv_tb_raw, merge_config, snp_vcf_gr) %>%
        expect_equal(as_granges(merged_tb))
    
    # test relative distance only
    merge_config <- list(
        # no limit
        maximum.gap.allowed = NA,
        call.extension.percent = 30,
        merge.gap.absolute = 0,
        merge.gap.snps = 0
    )
    df.or.GR <- tibble(
        seqnames = 'chr1',
        start = c(1000, 8000, 21000, 28000, 36000, 50000, 60000, 1e7, 1.6e7),
        width = c(5000, 5000,  5000,  5000,  5000,  5000,  5000, 5e6, 1e5),
        sample_id = 'test', 
        CNV_caller = 'test',
        CNV_type = 'test',
        CN = 1,
        ID = paste0('test', 1:9)
    )
    # empty_snps <- GRanges(sample_id = character(), FILTER = character(), LRR = numeric())
    expected <- tibble(
        seqnames = 'chr1',
        # 1+1, 1+1+1, 1, 1 1+1
        start = c( 1000, 21000, 50000, 60000, 1e7) %>% as.integer(),
        width = c(12000, 20000, 5000, 5000, 6.1e6) %>% as.integer(), 
        sample_id = 'test', 
        CNV_caller = 'test',
        CNV_type = 'test',
        n_initial_calls = c(2, 3, 1, 1, 2),
        initial_call_details = c('NA_CN1;NA_CN1', 'NA_CN1;NA_CN1;NA_CN1', NA_character_, NA_character_, 'NA_CN1;NA_CN1'),
        CN = 1,
    ) %>% as_granges() %>%
        mutate(
            ID = paste(CNV_caller, CNV_type, seqnames, start, end, sep='_'),
            n_probes = c(34,4,1,1,2),
            n_uniq_probes = c(33,4,1,1,2),
            probe_density_Mb = n_uniq_probes / width * 1e6,
        )
    merge_calls(df.or.GR, merge_config, snp_vcf_gr) %>%
        expect_equal(expected)
    # reduce % & test that max.gap is applied correctly
    merge_config <- list(
        maximum.gap.allowed = 10000,
        call.extension.percent = 20,
        merge.gap.absolute = 0,
        merge.gap.snps = 0
    )
    expected <- tibble(
        seqnames = 'chr1',
        # 1+1, 1+1, 1, 1, 1 1, 1
        start = c( 1000, 21000, 36000, 50000, 60000, 1e7, 1.6e7),
        width = c(12000, 12000,  5000,  5000,  5000, 5e6, 1e5),
        sample_id = 'test', 
        CNV_caller = 'test',
        CNV_type = 'test',
        n_initial_calls = c(2, 2, 1, 1, 1, 1, 1),
        initial_call_details = c('NA_CN1;NA_CN1', 'NA_CN1;NA_CN1', NA, NA, NA, NA, NA),
        CN = 1,
    ) %>% as_granges() %>%
        mutate(
            ID = paste(CNV_caller, CNV_type, seqnames, start, end, sep='_'),
            n_probes = c(34,3,1,1,1,1,1),
            n_uniq_probes = c(33,3,1,1,1,1,1),
            probe_density_Mb = n_uniq_probes / width * 1e6,
        )
    merge_calls(df.or.GR, merge_config, snp_vcf_gr) %>%
        expect_equal(expected)
      
    # Test for SNP based call merging
    merge_config <- list(
        maximum.gap.allowed = NA,
        call.extension.percent = 0,
        merge.gap.absolute = 0,
        merge.gap.snps = 5
    )
    df.or.GR <- tibble(
        seqnames = 'chr1',
        start = c( 100, 3000, 7400, 9600, 1e6, 2e6-1, 2875400),
        end = c(c(200, 4000, 8400, 10000)-1, 1250000, 2000400, 3000399),
        #width = c( 100, 1000, 1000,  400, 125e4-1, 401, 125e4),
        sample_id = 'test', 
        CNV_caller = 'test',
        CNV_type = 'test',
        CN = 1,
        ID = paste0('test', 1:7)
    )
    expected <- tibble(
        seqnames = 'chr1',
        # 1+1, 1+1, 1, 1, 1
        start = c( 100, 7400, 1e6, 2e6-1, 2875400),
        #width = c(3900, 2600, 125e4-1, 401, 125e4),
        end = c(c(4000, 10000)-1, 1250000, 2000400, 3000399),
        sample_id = 'test', 
        CNV_caller = 'test',
        CNV_type = 'test',
        n_initial_calls = c(2, 2, 1, 1, 1),
        initial_call_details = c('NA_CN1;NA_CN1', 'NA_CN1;NA_CN1', NA, NA, NA),
        CN = 1,
    ) %>% as_granges() %>%
        mutate(
            ID = paste(CNV_caller, CNV_type, seqnames, start, end, sep='_'),
            n_probes = c(14, 15, 3, 2, 3),
            n_uniq_probes = c(13, 15, 3, 2, 3),
            probe_density_Mb = n_uniq_probes / width * 1e6
        )
    merge_calls(df.or.GR, merge_config, snp_vcf_gr) %>%
        expect_equal(expected)
    # test with larger allwoed SNP gap (6 instead of 5)
    merge_config$merge.gap.snps <- 6
    expected <- tibble(
        seqnames = 'chr1',
        # 1+1, 1+1, 1+1+1
        start = c( 100, 7400, 1e6),
        #width = c(3900, 2600, 2e6+400),
        end = c(c(4000, 10000)-1, 3000399),
        sample_id = 'test', 
        CNV_caller = 'test',
        CNV_type = 'test',
        n_initial_calls = c(2, 2, 3),
        initial_call_details = c('NA_CN1;NA_CN1', 'NA_CN1;NA_CN1', 'NA_CN1;NA_CN1;NA_CN1'),
        CN = 1,
    ) %>% as_granges() %>%
        mutate(
            ID = paste(CNV_caller, CNV_type, seqnames, start, end, sep='_'),
            n_probes = c(14, 15, 20),
            n_uniq_probes = c(13, 15, 20),
            probe_density_Mb = n_uniq_probes / width * 1e6
        )
    merge_calls(df.or.GR, merge_config, snp_vcf_gr) %>%
        expect_equal(expected)
    
    # Test combination of modes, incl max.gap
    merge_config <- list(
        maximum.gap.allowed = 10000,
        call.extension.percent = 20,
        merge.gap.absolute = 1000,
        merge.gap.snps = 4
    )
    df.or.GR <- tibble(
        seqnames = 'chr1',
        # gains
        # 1+2: snp [3], 2+3: rel.distance [2000]
        # 4+5: !snps[5], rel [10k]  & =max.gap [10k]
        # 6,7: snps[1]|rel[200k] & >max.gap
        # losses
        # 8+9: abs [950]
        # 10+11: snps[4]|rel[10k-1]  & <max.gap [10k-1]
        # LOH
        # 12,13,14: none of the criteria
        start = c(
            c(100, 3000, 6000, 8e4, 14e4, 2e6+400, 2750400),
            c(  1, 1000, 8e4, 14e4), c(1, 3999, 21999)
        ),
        width = c(
            c(100, 1000, 16000, 5e4,  5e4, 5e5+1, 6e5+1),
            c( 49,  600, 5e4+1, 5e4), c(1000, 3001, 5e4+2)
        ),
        sample_id = 'test',
        CNV_caller = 'test',
        CNV_type = c(rep('DUP', 7), rep('DEL', 4), rep('LOH', 3)),
        CN = c(rep(3, 7), rep(1, 4), rep(0, 3)),
        ID = paste0('test', 1:14)
    )
    expected <- tibble(
        seqnames = 'chr1',
        start = c(  100,  8e4, 2e6+400, 2750400,    1,  8e4,    1, 3999, 21999),
        width = c(21900, 11e4,   5e5+1,   6e5+1, 1599, 11e4, 1000, 3001, 5e4+2),
        sample_id = 'test', 
        CNV_caller = 'test',
        CNV_type = c('DUP', 'DUP', 'DUP', 'DUP', 'DEL', 'DEL', 'LOH', 'LOH', 'LOH'),
        n_initial_calls = c(3, 2, 1, 1, 2, 2, 1,1,1),
        initial_call_details = c('NA_CN3;NA_CN3;NA_CN3', 'NA_CN3;NA_CN3', NA, NA, 'NA_CN1;NA_CN1', 'NA_CN1;NA_CN1', NA, NA, NA),
        CN = c(rep(3, 4), rep(1, 2), rep(0, 3)),
    ) %>% as_granges() %>%
        mutate(
            ID = paste(CNV_caller, CNV_type, seqnames, start, end, sep='_'),
            n_probes = c(40, 18, 5, 5, 10, 18, 8, 11, 7),
            n_uniq_probes = c(39, 18, 5, 5, 10, 18, 8, 11, 7),
            probe_density_Mb = n_uniq_probes / width * 1e6
        )
    merge_calls(df.or.GR, merge_config, snp_vcf_gr) %>%
        expect_equal(expected)
})

test_that("merge_calls & add_call_prefilters", {
    tool_config <- list(
        filter.minlength = 1000,
        filter.minprobes = 5,
        filter.mindensity.Mb = 10
    )
    merge_config <- list(
        maximum.gap.allowed = 1000,
        call.extension.percent = 0,
        merge.gap.absolute = 500,
        merge.gap.snps = 0
    )
    merge_calls(cnv_tb_raw, merge_config, snp_vcf_gr) %>%
        add_call_prefilters(tool_config) %>%
        expect_equal(as_granges(merged_filtered_tb))  
} )


test_that("test empty calls", {
    tool_config <- list(
        filter.minlength = 1000,
        filter.minprobes = 5,
        filter.mindensity.Mb = 10
    )
    merge_config <- list(
        maximum.gap.allowed = 1000,
        call.extension.percent = 0,
        merge.gap.absolute = 500,
        merge.gap.snps = 0
    )
    
    empty_raw <- cnv_tb_raw %>%
        filter(sample_id == 'non_existent_sample')
    empty_snps <- cnv_tb_snps %>%
        filter(sample_id == 'non_existent_sample')
    empty_merged <- merged_tb %>%
        filter(sample_id == 'non_existent_sample') %>%
        as_granges()
    
    as_granges(empty_raw) %>%
        add_snp_probe_counts(snp_vcf_gr) %>%
        expect_equal(as_granges(empty_snps))
    
    merge_calls(empty_raw, merge_config, snp_vcf_gr) %>%
        expect_equal(empty_merged)
    
    empty_merged$FILTER <- character()
    merge_calls(empty_raw, merge_config, snp_vcf_gr) %>%
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
        call.merging = list(
            maximum.gap.allowed = 1000,
            call.extension.percent = 0,
            merge.gap.absolute = 500,
            merge.gap.snps = 0
        ),
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
                snp_vcf_gr, 
                tool_config
            )
    )
    
    snp_vcf_gr <- fix_CHROM_format(snp_vcf_gr, 'NCBI')
    as_granges(cnv_tb_raw) %>%
            fix_CHROM_format('NCBI') %>%
            apply_preprocessing(
                snp_vcf_gr, 
                tool_config
            ) %>% 
        expect_equal(expected)
})

test_that('get_median_LRR', {
    
    expected_LRR <- merged_tb %>%
        mutate(LRR = c(1, 1.3, 0.89, 2, 1.385, 0, -0.88, -1.31, -0.895))
    
    merged_tb %>%
        as_granges() %>%
        get_median_LRR(snp_vcf_gr) %>%
        expect_equal(as_granges(expected_LRR))
    
})

# test BAF cluster determination
# Need to add cases where this should fail
#BAF clusters:
#4, 2, 4, 3, 4, 2, 2, 2, 2  