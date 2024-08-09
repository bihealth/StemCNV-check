library(testthat)

library(tidyverse)
library(plyranges)
library(vcfR)

source(test_path('../../stemcnv_check/scripts/R/vcf_io_functions.R'))

#  Functions to test:
# - vcfR_to_tibble
# - get_fix_section
# - get_gt_section
# - parse_cnv_vcf
# Note: static_cnv_vcf_header doesn't really need a test

snp_vcf <- read.vcfR(test_path('../data/minimal_probes.vcf'))

test_that("vcfR_to_tibble", {    
    expected_tb <- tibble(
        CHROM = 'chr1',
        POS = c(100, 105, 110, 115, 115, 199, 1000, 1300, 1599, 3000, 3250, 3500, 3500, 3750, 3999) %>% as.integer(),
        ID = c(NA, NA, NA, NA, 'dummy', NA, NA, NA, NA, NA, NA, NA, 'DUP', NA, NA),
        REF = NA_character_, # or "."
        ALT = NA_character_, # or "."
        QUAL = NA_real_, # or "."
        FILTER = c(rep('PASS', 4), 'LOWQUAL', rep('PASS', 10)),
        GenTrain_Score = c(rep(0.8, 4), 0.2, rep(0.8, 7), 0.9, 0.8, NA),
        CSQ = c('gene_id|gene_name', NA, NA, NA, NA, NA, NA, NA, NA, 'ABC|abc', NA, NA, NA, NA, NA),
        sample_id = 'test_sample',
        GT = NA_character_, # or "."
        LRR = c(1, 1.1, 0.7, 1, 1, 0.9, 1.3, 1.2, 1.3, 0.85, 0.82, 0.8, 1.15, 0.87, 0.88),
        BAF = c(1, 0, 0.3, 0.7, 1, 1, 1, 1, 0, 0.33, 1, 1, 0.67, 1, 0.67)        
    )
    
    vcfR_to_tibble(snp_vcf) %>%
        filter(POS < 4000) %>% 
        expect_equal(expected_tb)
})

# test_that('static_cnv_vcf_header', {
#     tool_config <- list(
#         filter.minsize = 1000,
#         filter.minprobes = 5,
#         filter.mindensity = 10
#     )
# })

#TODO: add additional tb with fully processed/annotated CNV calls
# > should probablt move those fixtures to testthat.R
cnv_tb <- tibble(
  seqnames  = rep('chr1', 9),
  start = c(100, 1000, 3000, 6000, 1.0e8, 7400,  9000, 12e4, 1e6) %>% as.integer(),
  end   = c(200, 1600, 5400, 7000, 1.1e8, 8400, 10000, 15e4, 3e6+400)-1 %>% as.integer(),
  width = c(100, 600, 2400, 1000, 1e7, 1000, 1000, 3e4, 2e6+400),
  sample_id = 'test_sample',
  CNV_caller = 'Test',
  CNV_type = c(rep('DUP', 5), rep('DEL', 4)),
  n_initial_calls = c(1, 1, 2, 1, 1, 1, 2, 3, 2),
  initial_call_details = c(NA, NA, '3000_3999_CN3,4400_5399_CN3', NA, NA, 
                           NA, '9000_9399_CN1,9600_9999_CN1', '120000_129999_CN1,130000_139999_CN0,140000_149999_CN1',
                            '1000000_1999999_CN1,2000400_3000399_CN1'),
  CN = c(rep(3, 5), rep(1, 4)),  
  ID = paste(CNV_caller, CNV_type, seqnames, start, end, sep='_'),
  n_probes =      c(5, 3, 11, 5, 20, 5, 10, 15, 20),
  n_uniq_probes = c(5, 3, 10, 5, 20, 5, 10, 15, 20),
  probe_density_Mb = n_uniq_probes / width * 1e6,
  FILTER = c('Size', 'Size;n_probes', 'PASS', 'PASS', 'Density', 'PASS', 'PASS', 'PASS', 'Density')
)

test_that('get_fix_section', {
    expected_fix <- tibble(
        CHROM = 'chr1',
        POS = cnv_tb$start - 1,
        ID = cnv_tb$ID,
        REF = '.',
        ALT = paste0('<', cnv_tb$CNV_type, '>'),
        QUAL = '.',
        FILTER = cnv_tb$FILTER,
        # 'END={end};SVLEN={width};SVCLAIM=D;N_PROBES={n_probes};N_UNIQ_PROBES={n_uniq_probes};PROBE_DENS={probe_density_Mb}'
        INFO = paste0(
            str_glue('END={cnv_tb$end};SVLEN={cnv_tb$width};'),
            'SVCLAIM=D;',
            str_glue('N_PROBES={cnv_tb$n_probes};N_UNIQ_PROBES={cnv_tb$n_uniq_probes};'),
            str_glue('PROBE_DENS={round(cnv_tb$probe_density_Mb, 2)}')            
        ),
    ) %>%
        as.matrix()
    
    get_fix_section(cnv_tb) %>%
        expect_equal(expected_fix)

})

#TODO: test different CNs (LOH, 0, >=4,, X&Y on male)
test_that('get_gt_section', {
    # FORMAT keys: GT, CN, TOOL (str desc), LRR (median)
    expected_gt <- tibble(
        FORMAT = 'GT:CN:TOOL:LRR',
        test_sample = paste(
            '0/1',
            cnv_tb$CN,
            paste0(
                "caller=Test;",
                "n_initial_calls=", cnv_tb$n_initial_calls, ";",
                "initial_call_details=", ifelse(is.na(cnv_tb$initial_call_details), '.', cnv_tb$initial_call_details)
            ),
            #expcted LRR medians from minimal_porbes vcf
            c(1, 1.3, 0.89, 2, 1.385, 0, -0.88, -1.31, -0.895),
            #BAF clusters:
            #4, 2, 4, 3, 4, 2, 2, 2, 2
            sep = ":"
        )
    ) %>% as.matrix()
    
    get_gt_section(cnv_tb, snp_vcf) %>%
        expect_equal(expected_gt)

})