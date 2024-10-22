library(testthat)

library(tidyverse)
library(plyranges)

source(test_path('../../stemcnv_check/scripts/R/helper_functions.R'))
source(test_path('../../stemcnv_check/scripts/R/processCNVs_annotate_check-score.R'))

# Functions:
# - annotate_cnv.check.score
# - annotate_precision.estimates
# - annotate_Call.label

config <- list(
    'genome_version' = 'hg19',
    'global_settings' = list(
        'hg19_gtf_file' = test_path('../data/hg_minimal.gtf'),
        'hg19_genomeInfo_file' = test_path('../data/gr_info_minimal.tsv')
    ),
    'settings' = list(
        'CNV_processing' = list(
            'gene_overlap' = list(
                'exclude_gene_type_regex' = c(),
                'include_only_these_gene_types' = c()
            )
        ),
        'chromosomes' = paste0('chr', c(1:22, 'X', 'Y'))
    )
)
gtf_file <- test_path('../data/hg_minimal.gtf')
ginfo_file <- test_path('../data/gr_info_minimal.tsv')

base_tb <- tibble(
    seqnames = 'chr1',
    start = c(4000, 10000, 40000, 28000000, 28060000, 40000, 5000000, 3000) %>% as.integer(),
    end   = c(5500, 14000, 50000, 28055000, 28065000, 50000, 7000000, 60000) %>% as.integer(),
    sample_id = 'test_sample',
    CNV_type = c('gain', 'gain', 'gain', 'gain', 'loss', 'loss', 'LOH', 'LOH'),
    ID = paste('combined', CNV_type, seqnames, start, end, sep='_'),
    CNV_caller = c('StemCNV-check', 'toolA', 'toolA', 'StemCNV-check', 'StemCNV-check', 'toolA', 'toolB', 'toolB'),
    n_probes = c(15, 100, 100, 150, 100, 100, 50, 50),
    n_uniq_probes = c(15, 100, 100, 150, 100, 100, 50, 50),
    CN = c(3, 3, 3, 4, 1, 1, 2, 2),
    FILTER = c('min_size', 'Probe_dens;probe_gap', NA_character_, 'high_probe_dens', 'test-dummy;high_probe_dens', 'probe_gap', NA_character_, NA_character_),
    reference_overlap = c(T, T, F, F, F, F, F, T),
    reference_coverage = c(100, 85.01, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, 60),
    reference_caller = c('StemCNV-check', 'faketool', NA_character_, NA_character_, NA_character_,NA_character_,NA_character_, 'toolA'),
    stemcell_hotspot = c(NA, NA, 'chr1:40000-50000', 'dummyC', NA, '1p36|chr1:40000-50000', NA, NA),
    dosage_sensitive_gene = c(NA, NA, 'dummyB', NA, 'dummyC', NA, NA, NA),
    cancer_gene = c(NA, 'DDX11L1', NA, NA, NA, NA, NA, 'DDX11L1'),
    ROI_hits = c('fake-ROI', NA, NA, NA, 'dummyC', NA, NA, NA),    
    Gap_percent = c(0, 2000/4001, 2000/10001, 25000/55001, 1000/5001, 2000/10001, 1e6/(2e6+1), 0),
    probe_coverage_gap = c(FALSE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE),
    high_probe_density = c(NA, NA, NA, TRUE, TRUE, NA, FALSE, FALSE)
)

expected_gene_tb <- base_tb %>%
    mutate(
        width = c( 1501, 4001, 10001, 55001, 5001, 10001, 2000001, 57001) %>% as.integer(),
        n_genes = c(1, 1, 1, 1, 1, 1, 0, 3)  %>% as.integer(),
        overlapping_genes = c('dummyA', 'DDX11L1', 'dummyB', "dummyC", "dummyC", 'dummyB', NA, 'dummyA|DDX11L1|dummyB'),
    )

expected_final_tb <- expected_gene_tb %>%
    mutate(
        seqnames = factor(seqnames, levels = genomeStyles('Homo_sapiens')[['UCSC']])
    )

test_that("annotate_gene_overlap", {
    # Test scenarios:
    # - CNV has no gene overlap
    # - CNV has (partial) gene overlap
    # - CNV has multiple genes overlapping    
    gr_genes <- load_gtf_data(gtf_file, config)
    
    annotate_gene_overlaps(as_granges(base_tb), gr_genes) %>%
        expect_equal(as_granges(expected_gene_tb))
    
})


test_that("Annotate CNV check scores", {
    CNV_size_score <- function(len) {   1/3 * log(len) * log(len) - 15 }
    LOH_size_score <- function(len) { 0.275 * log(len) * log(len) - 15 }
    config$settings$CNV_processing$Check_score_values <- list(
        'any_roi_hit' = 50,
        'any_other_gene' = 0.2,
        'large_CN_size_modifier' = 1.5
    )
  
    stemcell_hotspot_gr <- tibble(
        seqnames = 'chr1',
        start = c(28050000, 40000, 0),
        end = c(28070000, 50000, 7200000),
        strand = c('+', '*', '*'),
        list_name = 'test-list',
        hotspot = c('dummyC', 'chr1:40000-50000', '1p36'),
        mapping = c('gene_name', 'position', 'gband'),
        call_type = c('gain', 'any', 'loss'),
        check_score = c(15, 30, 10),
        source = c('dummy', 'dummy', 'dummy'),
        comment = NA
    ) %>% as_granges()
  
    cancer_gene_gr <- tibble(
        seqnames = 'chr1',
        start = 11873,
        end = 14409,
        strand = '+',
        list_name = 'test-list',
        hotspot = 'DDX11L1', 
        mapping = 'gene_name', 
        call_type = 'any', 
        check_score = 5,
        source = 'dummy',
        comment = NA
    ) %>% as_granges()
    
    dosage_sensitive_gene_gr <- tibble(
        seqnames = 'chr1',
        start = c(4000, 20000, 28050000),
        end = c(5000, 50000, 28070000),
        strand = '+',
        list_name = 'test-dosage',
        hotspot = c('dummyA', 'dummyB', 'dummyC'),
        mapping = 'gene_name', 
        call_type = c('loss', 'gain', 'loss'),
        check_score = 7,
        source = 'dummy',
        comment = NA
    ) %>% as_granges()
  
    expected_tb <- expected_final_tb %>%
        mutate(
            # Test scenarios:
            Check_Score = c(
                # ROI (50) + 1 other gene w/o hotspot 
                CNV_size_score(1501) + 50 + 0.2, 
                # cancer gene (5)
                CNV_size_score(4001) + 5, 
                # hotspot + dosage gene (30 & 7)
                CNV_size_score(10001) + 30 + 7,
                # CN4, hotspot gene (15)
                CNV_size_score(55001) * 1.5 + 15, 
                # ROI hit (50) + dosage gene (7)
                CNV_size_score(5001) + 50 + 7,
                # 2 hotspots (30 & 10) + 1 other gene
                CNV_size_score(10001) + 30 + 10 + 0.2,
                # no genes
                LOH_size_score(2000001) + 0,
                # cancer gene (5) + 2 other genes
                LOH_size_score(57001) + 5 + 0.4
            )
        )
    
    annotate_cnv.check.score(
        expected_final_tb,
        stemcell_hotspot_gr,
        dosage_sensitive_gene_gr,
        cancer_gene_gr,
        config$settings$CNV_processing$Check_score_values
    ) %>%
        expect_equal(expected_tb)
})



# test_that("Annotate precision estimates", {
# - CBS only
# - PennCNV only
# - CBS + PennCNV
# - different sizes


test_that("Annotate call label", {
    # Test scenarios:
    # - ref GT ( = ref coverage >= X% ?!) [1,2,8]
    # - critical score [5]
    # - reportable score [4]
    # - critical score, but crit. excl list (-> reportable)  [6]
    # - reportable score, but excl list (-> NA) [added, 9]
    # - NA [7]
    call_cat_config <- list(
        check_score.critical = 53,
        filters.exclude.critical = c('probe_gap'),
        check_score.reportable = 50,
        filters.exclude.reportable = c('test-dummy')
    )  
    
    input_tb <- expected_final_tb %>%
        mutate(
            Check_Score = c(53.03098, 12.93180, 50.27740, 52.06978, 66.18200, 53.47740, 42.88782, 23.37815)
        ) %>%
        bind_rows(
           expected_final_tb[6,] %>% mutate(FILTER = 'test-dummy')
        )
    
    expected_tb <- input_tb %>%
        mutate(
            Call_label = c('Reference genotype', 'Reference genotype', 'Reportable', 'Reportable', 'Critical', 'Reportable', NA, 'Reference genotype', NA)
        ) 
    expect_equal(annotate_call.label(input_tb, call_cat_config), expected_tb)
    
    # No reportable defined, only 'downgraded' critical becomes reportable
    call_cat_config <- list(
        check_score.critical = 55,
        filters.exclude.critical = c('probe_gap'),
        check_score.reportable = NULL,
        filters.exclude.reportable = c()
    )  
    
    expected_tb <- input_tb %>%
        mutate(
            Call_label = c('Reference genotype', 'Reference genotype', NA, NA, 'Critical', NA, NA, 'Reference genotype', NA)
        ) 
    expect_equal(annotate_call.label(input_tb, call_cat_config), expected_tb)
})
