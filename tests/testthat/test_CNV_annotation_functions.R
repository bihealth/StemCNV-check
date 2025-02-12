library(testthat)

library(tidyverse)
library(plyranges)

source(test_path('../../stemcnv_check/scripts/R/helper_functions.R'))
source(test_path('../../stemcnv_check/scripts/R/hotspot_functions.R'))
source(test_path('../../stemcnv_check/scripts/R/CNV_annotation_functions.R'))


config <- list(
    'snakedir' = '',
    'genome_version' = 'hg19',
    'global_settings' = list(
        'hg19_gtf_file' = test_path('../data/hg_minimal.gtf'),
        'hg19_genomeInfo_file' = test_path('../data/gr_info_minimal.tsv')
    ),
    'settings' = list(
        'CNV_processing' = list(
            'gene_overlap' = list(
                'exclude_gene_type_regex' = c(),
                'include_only_these_gene_types' = c(), # c('lncRNA', 'miRNA', 'protein_coding'),
                'whitelist_hotspot_genes' = FALSE,
                'stemcell_hotspot_list' = test_path('../data/minimal-hotspots.tsv')
            )
        ),
        'vcf_output' = list(
            'chrom_style' = 'keep-original'
        )
    )
)

gr_info <- load_genomeInfo(config$global_settings$hg19_genomeInfo_file, config)
gr_genes <- load_gtf_data(config$global_settings$hg19_gtf_file, config)

sample_cnvs <- tibble(
  seqnames = 'chr1',
  start = c(4000, 10000, 28000000, 28060000, 40000, 5000000, 31000000),
  end   = c(5500, 14000, 28055000, 28065000, 50000, 7000000, 32000000),
  sample_id = 'test_sample',
  CNV_type = c('gain', 'gain', 'gain', 'loss', 'loss', 'loss', 'loss'),
  ID = paste('combined', CNV_type, seqnames, start, end, sep='_'),
  CNV_caller = c('StemCNV-check', 'toolA', 'StemCNV-check', 'StemCNV-check', 'toolA', 'toolB', 'toolB'),
  n_probes = c(15, 100, 150, 250, 100, 50, 50),
  n_uniq_probes = c(15, 100, 150, 100, 100, 50, 50),
  CN = c(3, 3, 4, 1, 1, 1, 0),
  FILTER = c('min_size', 'min_density', NA_character_, 'test-dummy', NA_character_, NA_character_, NA_character_),
  reference_overlap = c(T, T, F, F, F, F, T),
  reference_coverage = c(100, 85.01, NA_real_, NA_real_, NA_real_, NA_real_, 60) %>% as.list(),
  reference_caller = c('StemCNV-check', 'faketool', NA_character_, NA_character_,NA_character_,NA_character_, 'toolA'),
  test_hits = c(NA, 'DDX11L1', 'dummyC', NA, '1p36,chr1:40000-50000', '1p36', '1p35.2')
) %>% as_granges()

### Impact Lists ###

# 1 - not hit
# 2 - gene hit, any (but not gband): DDX11L1
# 3 - gene hit, gain: dummyC
# 4 - gene hit, no match
# 5 - position hit (any) + gband_hit (loss): chr1:40000-50000 + 1p36
# 6 - gband hit (any subband): 1p36
# 7 - gband hit (specific subband): 1p35.2

#annotate_impact_lists
test_that('annotate_impact_lists', {
    
    hotspots <- parse_hotspot_table(read_tsv(test_path('../data/minimal-hotspots.tsv')), gr_genes, gr_info)
    
    expected_gr <- sample_cnvs %>%
        mutate(test = c(NA, 'DDX11L1', 'dummyC', NA, '1p36|chr1:40000-50000', '1p36', '1p35.2'))
    expect_equal(annotate_impact_lists(sample_cnvs, hotspots, 'test'), expected_gr)

    # test empty hotspots
    expected_gr <- sample_cnvs %>%
        mutate(test = NA_character_)
    expect_equal(annotate_impact_lists(sample_cnvs, GRanges(), 'test'), expected_gr)

})


test_that('annotate_roi', {
    
    sample_tb <- tibble(
        Sample_ID = 'test_sample',
        Regions_of_Interest = 'chr1:10000-20000;nametest|chr1:5000000-5500000;chrX:30000-40000;1p35'
    )

    expected_gr <- sample_cnvs %>%
        mutate(ROI_hits = c(NA, 'chr1:10000-20000', '1p35', '1p35', NA, 'chr1:5000000-5500000', '1p35'))
    roi_tb <- get_roi_tb('test_sample', sample_tb, config)
    
    expect_equal(annotate_roi(sample_cnvs, roi_tb, gr_genes, gr_info, config), expected_gr)
    
    # test empty roi definition
    expected_gr <- sample_cnvs %>%
        mutate(ROI_hits = NA_character_)
    sample_tb$Regions_of_Interest <- ''
    roi_tb <- get_roi_tb('test_sample', sample_tb, config)
    expect_equal(annotate_roi(sample_cnvs, roi_tb, gr_genes, gr_info, config), expected_gr)
})


### Array Features ###

# Gaps:
# 1 - no gaps in CNV
# 2 - single gap in CNV         -> above th, true
# 3 - multiple gaps in one CNV  -> above th, true
# 4 - CNV has single gap,    but below min.perc.gap_area
# 5 - CNV has multiple gaps, but below min.perc.gap_area
# 6 - CNV has gaps, but ! (gap_slope * Gap_percent + gap_intercept) <= log2(n_uniq_probes)
# 7 - no gaps in CNV
# Extra tests:
# - [x] gap overlaps CNV border -> error
# - [ ] gap fully conatins CNV  -> error

test_that("Annotate CNVs with gaps", {
  gapfile <- test_path('../data/gaps_minimal.bed')
  gap_area.uniq_probes.rel <- list(-12, 12.5)
  min.perc.gap_area <- 0.33
  
  expexted_gr <- sample_cnvs %>%
    mutate(
      Gap_percent = c(0, 2000/4001, 25000/55001, 1000/5001, 2000/10001, 1e6/(2e6+1), 0),
      probe_coverage_gap = c(FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE),
      FILTER = c('min_size', 'min_density;probe_gap', 'probe_gap', 'test-dummy', NA_character_, NA_character_, NA_character_)
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

# Density:
# 1 - 
# 2 - 
# 3 -
# 4 - 
# 5 -
# 6 - 
# 7 - 

test_that("Annotate CNVs with probe density flags", {
  density_file <- test_path('../data/density_minimal.bed')
  density.quantile.cutoff <- 0.99
  
  expexted_gr <- sample_cnvs %>%
    mutate(
      high_probe_density = c(NA, NA, TRUE, TRUE, NA, FALSE, FALSE),
      FILTER = c('min_size', 'min_density', 'high_probe_dens', 'test-dummy;high_probe_dens', NA_character_, NA_character_, NA_character_)
    )
  
  annotate_high_density(sample_cnvs, density_file, density.quantile.cutoff) %>%
    expect_equal(expexted_gr)        

} )

### Genes & Check Score ###


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
    
    annotate_gene_overlaps(as_granges(base_tb), gr_genes) %>%
        expect_equal(as_granges(expected_gene_tb))
    
})


test_that("Annotate CNV check scores", {
    CNV_double_copy_score <- function(len) { 0.5 * log(len) * log(len) - 15 }
    CNV_single_copy_score <- function(len) { 0.333 * log(len) * log(len) - 15 }
    LOH_size_score <- function(len) { 0.275 * log(len) * log(len) - 15 }
    config$settings$CNV_processing$Check_score_values <- list(
        'any_roi_hit' = 50,
        'any_other_gene' = 0.2,
        'single_copy_factor' = 0.333,
        'double_copy_factor' = 0.5,
        'neutral_copy_factor' = 0.275,
        'flat_decrease' = 15
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
                CNV_single_copy_score(1501) + 50 + 0.2, 
                # cancer gene (5)
                CNV_single_copy_score(4001) + 5, 
                # hotspot + dosage gene (30 & 7)
                CNV_single_copy_score(10001) + 30 + 7,
                # CN4, hotspot gene (15)
                CNV_double_copy_score(55001) + 15, 
                # ROI hit (50) + dosage gene (7)
                CNV_single_copy_score(5001) + 50 + 7,
                # 2 hotspots (30 & 10) + 1 other gene
                CNV_single_copy_score(10001) + 30 + 10 + 0.2,
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
        config$settings$CNV_processing$Check_score_values,
        "f"
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
    # - reportable score, but excl list (-> basic) [added, 9]
    # - basic [7]
    call_cat_config <- list(
        Critical = list(
            minimum_check_score = 53,
            not_allowed_vcf_filters = list('probe_gap'),
            reference_match = FALSE
        ),
        Reportable = list(
            minimum_check_score = 50,
            not_allowed_vcf_filters = list('test-dummy'),
            reference_match = FALSE
        ),
        `de-novo` = list(
            minimum_check_score = 0,
            not_allowed_vcf_filters = list(),
            reference_match = FALSE
        ),
        `Reference genotype` = list(
            minimum_check_score = 0,
            not_allowed_vcf_filters = list(),
            reference_match = TRUE
        )
    )
    
    input_tb <- expected_final_tb %>%
        mutate(
            Check_Score = c(53.03098, 12.93180, 50.27740, 52.06978, 66.18200, 53.47740, 42.88782, 23.37815)
        ) %>%
        bind_rows(
           expected_final_tb[6,] %>% mutate(FILTER = 'test-dummy', Check_Score = 49.9)
        )
    
    expected_tb <- input_tb %>%
        mutate(
            Call_label = c(
                'Reference genotype', 'Reference genotype', 'Reportable', 'Reportable', 'Critical', 
                'Reportable', 'de-novo', 'Reference genotype', 'de-novo'
            )
        ) 
    expect_equal(annotate_call.label(input_tb, call_cat_config), expected_tb)
    
    # Same Check_Score, but different vcf filters allowed for Critical & Reportable
    # > values between 50 and 53 should be assigned to basic now
    call_cat_config$Reportable$minimum_check_score <- 53
    
    expected_tb <- input_tb %>%
        mutate(
            Call_label = c(
                'Reference genotype', 'Reference genotype', 'de-novo', 'de-novo', 'Critical',
                'Reportable', 'de-novo', 'Reference genotype', 'de-novo'
            )
        ) 
    expect_equal(annotate_call.label(input_tb, call_cat_config), expected_tb)
})

