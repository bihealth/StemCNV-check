library(testthat)

library(tidyverse)
library(plyranges)
library(yaml)


source(test_path('../../stemcnv_check/scripts/R/helper_functions.R'))
source(test_path('../../stemcnv_check/scripts/R/vcf_io_functions.R'))
source(test_path('../../stemcnv_check/scripts/R/CNV_preprocess_functions.R'))
source(test_path('../../stemcnv_check/scripts/R/CNV_comparison_functions.R'))

defined_labels <- yaml.load_file(test_path('../../stemcnv_check/control_files/label_name_definitions.yaml'))

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
    start = c(1000, 4000, 6000,  8000, 120000, 15e5, 18e5, 22.2e4, 27e4, 31e4, 36e4, 50e4, 55e4) %>% as.integer(),
    end   = c(2000, 5500, 7500, 10000, 140000, 17e5, 20e5, 23e4, 27.8e4, 35e4, 40e4, 54e4, 58e4) %>% as.integer(),
    CNV_type = c(rep('DUP', 3), rep('DEL', 10)),
    sample_id = 'test_sample',
    FILTER = c('min_size', 'Test;Test2', 'min_probes', rep(NA_character_, 10)),
    CNV_caller = 'toolA',
    n_probes = 10,
    n_uniq_probes = n_probes,
    probe_density_Mb = n_uniq_probes / (end - start + 1) * 1e6,
    n_initial_calls = 1,
    initial_call_details = NA_character_,
    CN = c(rep(3, 3), rep(1, 10)),
    ID = paste(CNV_caller, CNV_type, seqnames, start, end, sep='_'),
)

toolB <- tibble(
    seqnames = rep('chr1', 10),
    start = c(2000, 4000, 6000,  9700, 120000, 15e5, 21e4, 33e4, 52e4, 57e4) %>% as.integer(),
    end   = c(3000, 5500, 7500, 10000, 140000, 20e5, 30e4, 36e4, 56e4, 60e4) %>% as.integer(),
    CNV_type = c(rep('DUP', 2), rep('DEL', 8)),
    sample_id = 'test_sample',
    FILTER = c(NA_character_, NA_character_, 'min_density', rep(NA_character_, 7)),
    CNV_caller = 'toolB',
    n_probes = 5,
    n_uniq_probes = n_probes,
    probe_density_Mb = n_uniq_probes / (end - start + 1) * 1e6,
    n_initial_calls = 1,
    initial_call_details = NA_character_,
    CN = c(rep(3, 2), 1, 1, 0, rep(1, 5)),
    ID = paste(CNV_caller, CNV_type, seqnames, start, end, sep='_'),
)


combined_tools <- bind_ranges(
    #combined calls
    # sortebd by: CNV_type (DEL before DUP), start
    tibble(     
        seqnames = 'chr1',
        start = c(120000, 15e5, 4000) %>% as.integer(),
        end   = c(140000, 20e5, 5500) %>% as.integer(),
        sample_id = 'test_sample',
        CNV_type = c('DEL', 'DEL', 'DUP'),
        ID = paste(defined_labels$combined_cnvs, CNV_type, seqnames, start, end, sep='_'),
        n_initial_calls = c(2, 3, 2),
        initial_call_details = c(
            "toolA_120000-140000_CN1_cov100_PASS|toolB_120000-140000_CN0_cov100_PASS",
            "toolA_1500000-1700000_CN1_cov40_PASS|toolA_1800000-2000000_CN1_cov40_PASS|toolB_1500000-2000000_CN1_cov100_PASS",
            "toolA_4000-5500_CN3_cov100_Test&Test2|toolB_4000-5500_CN3_cov100_PASS"
        ),       
        CNV_caller = defined_labels$combined_cnvs,
        CN = c(0.5, 1, 3),
        # overlap_merged_call = NA_real_,
        # Recalculated based on vcf file; indiv call will keep whatever (fake) number they had before
        n_probes = c(11, 6, 5),
        n_uniq_probes = n_probes,
        # this *should* be correct
        probe_density_Mb = n_uniq_probes / (end - start + 1) * 1e6,
        FILTER = c(NA, NA, 'min_probes'),
    ) %>% as_granges(),
    # single calls
    # sorted by type, then caller, then start
    bind_rows(
        toolA[c(1,3,4,8:13),],
        toolB[c(1,3,4,7:10),]
    ) %>%
        arrange(CNV_type, CNV_caller, start) %>%
        mutate(
                n_initial_calls = 1,
                initial_call_details = NA_character_,
        ) %>%
        as_granges()
)

snp_vcf <- parse_snp_vcf(test_path('../data/minimal_probes.vcf'))

test_that("combine 2 CNV callers", {
    
    processing_config <- list(
        tool.overlap.greatest.call.min.perc = 50,
        tool.overlap.min.cov.sum.perc = 60,
        filter.minlength = 1000,
        # ensure one call gets the min_probes filer applied
        filter.minprobes = 6,
        filter.mindensity.Mb = 10
    )
  
    list(
        'toolA' = as_granges(toolA),
        'toolB' = as_granges(toolB)
    ) %>%
        bind_ranges() %>%
        combine_CNV_callers(processing_config, snp_vcf, defined_labels) %>% 
        expect_equal(combined_tools)        
} )

# add test with empty callset for one tool
test_that("No calls from 1 CNV caller", {
    processing_config <- list(
        tool.overlap.greatest.call.min.perc = 50,
        tool.overlap.min.cov.sum.perc = 60,
        filter.minlength = 1000,
        # ensure one call gets the min_probes filer applied
        filter.minprobes = 6,
        filter.mindensity.Mb = 10
    )
  
    list(
        'toolA' = as_granges(toolA),
        'toolB' = as_granges(toolB) %>% filter(sample_id == 'non_existent_sample')
    ) %>%
        bind_ranges() %>%
        combine_CNV_callers(processing_config, snp_vcf, defined_labels) %>%
        expect_equal(
            as_granges(toolA) %>%
                mutate(
                    initial_call_details = NA_character_,
                    n_initial_calls = 1
                )
        )        
})

# add test with completely empty callset
test_that("No calls at all", {
    processing_config <- list(
        tool.overlap.greatest.call.min.perc = 50,
        tool.overlap.min.cov.sum.perc = 60,
        filter.minlength = 1000,
        # ensure one call gets the min_probes filer applied
        filter.minprobes = 6,
        filter.mindensity.Mb = 10
    )
  
    empty_res <- list(
        'toolA' = as_granges(toolA) %>% filter(sample_id == 'non_existent_sample'),
        'toolB' = as_granges(toolB) %>% filter(sample_id == 'non_existent_sample')
    ) %>%
        bind_ranges() %>%
        combine_CNV_callers(processing_config, snp_vcf, defined_labels) 
    
    empty_gr <- as_granges(toolA) %>%
        arrange(CNV_type, CNV_caller, start) %>%
        mutate(
            initial_call_details = NA_character_,
            n_initial_calls = 1
        ) %>%
        filter(sample_id == 'non_existent_sample')
  
    #Somehow the *direct* empty Granges objects do not compare equal ??
    expect_equal(as_tibble(empty_res) %>% as_granges(), as_tibble(empty_gr) %>% as_granges())
} )

# Example calls for reference annotation testing

# Call with ref (100% match)
# Call with ref (85% match)
# Call with ref (wrong CNV type)
# Call with ref (ref <80% OV)
# Call with ref (call <80% OV)
# Call without ref
# Ref call without match

sample_cnvs <- tibble(   
  seqnames = 'chr1',
  start = c(4000, 10000, 12000, 15e4, 28e4, 35e4) %>% as.integer(),
  end   = c(5500, 11000, 14000, 20e4, 30e4, 40e4) %>% as.integer(),
  sample_id = 'test_sample',
  CNV_type = c('DUP', 'DUP', 'DEL', 'DEL', 'DEL', 'DEL'),
  ID = paste('combined', CNV_type, seqnames, start, end, sep='_'),
  CNV_caller = c(defined_labels$combined_cnvs, 'toolA', defined_labels$combined_cnvs, defined_labels$combined_cnvs, 'toolA', 'toolB'),
  n_probes = c(15, 5, 15, 25, 5, 5),
  CN = c(3, 3, 0, 1, 1, 1), 
)

ref_cnvs <- tibble(   
  seqnames = 'chr1',
  start = c(4000, 10000, 12000, 18e4, 25e4, 45e4) %>% as.integer(),
  end   = c(5500, 10850, 14000, 20e4, 30e4, 50e4) %>% as.integer(),
  sample_id = 'test_sample',
  CNV_type = c('DUP', 'DUP', 'DUP', 'DEL', 'DEL', 'DEL'),
  ID = paste('combined', CNV_type, seqnames, start, end, sep='_'),
  CNV_caller = c(defined_labels$combined_cnvs, 'faketool', defined_labels$combined_cnvs, defined_labels$combined_cnvs, 'toolA', 'toolB'),
  n_probes = c(15, 20, 15, 25, 5, 5),
  CN = c(3, 4, 3, 1, 1, 1),
) %>%
  as_granges()

test_that("annotate_reference_overlap", {
  min.reciprocal.coverage.with.ref <- 80
  
  expected_gr <- sample_cnvs %>%
    mutate(
        reference_overlap = c(T, T, F, F, F, F),
        reference_coverage = c(100, 85.01, NA_real_, NA_real_, NA_real_, NA_real_),
        reference_caller = c(defined_labels$combined_cnvs, 'faketool', NA_character_, NA_character_,NA_character_,NA_character_)
    ) %>%
    as_granges()
  
  annotate_reference_overlap(
      as_granges(sample_cnvs),
      ref_cnvs,
      min.reciprocal.coverage.with.ref
  ) %>% 
    expect_equal(expected_gr) 
} )

test_that("annotate_reference_overlap, empty ref", {
  min.reciprocal.coverage.with.ref <- 80
  
  expected_gr <- sample_cnvs %>%
    # This is somehow needed to successfully run the test via Rscript
    as_granges() %>% as_tibble() %>%
    mutate(
        reference_overlap = F, 
        reference_coverage = rep(NA_real_, 6),
        reference_caller = rep(NA_character_, 6)
    ) %>%
    as_granges()
  
  annotate_reference_overlap(
      as_granges(sample_cnvs), 
      ref_cnvs %>% filter(sample_id == 'non_existent_sample'),
      min.reciprocal.coverage.with.ref
  ) %>%
    expect_equal(expected_gr)        
} )

test_that("annotate_reference_overlap, empty input", {
  min.reciprocal.coverage.with.ref <- 80
  
  expected_gr <- sample_cnvs %>%
    # This is somehow needed to successfully run the test via Rscript
    as_granges() %>% as_tibble() %>%
    mutate(
        reference_overlap = F, 
        reference_coverage = rep(NA_real_, 6),
        reference_caller = rep(NA_character_, 6)
    ) %>%
    as_granges() %>%
    filter(sample_id == 'non_existent_sample')
  
  annotate_reference_overlap(
      sample_cnvs %>% as_granges() %>% filter(sample_id == 'non_existent_sample'),
      ref_cnvs, 
      min.reciprocal.coverage.with.ref
  ) %>%
    expect_equal(expected_gr)        
} )
