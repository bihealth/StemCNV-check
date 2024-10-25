library(testthat)

library(tidyverse)
library(plyranges)

source(test_path('../../stemcnv_check/scripts/R/helper_functions.R'))
source(test_path('../../stemcnv_check/scripts/R/hotspot_functions.R'))

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
                'include_only_these_gene_types' = c('lncRNA', 'miRNA', 'protein_coding'),
                'stemcell_hotspot_list' = test_path('../data/minimal-hotspots.tsv'),
                'cancer_gene_list' = test_path('../data/minimal-hotspots.tsv')
            ),
            'Check_score_values' = list('any_roi_hit' = 50)
        ),
        'vcf_output' = list(
            'chrom_style' = 'keep-original'
        )
    )
)
gr_info <- load_genomeInfo(config$global_settings$hg19_genomeInfo_file, config)
gr_genes <- load_gtf_data(config$global_settings$hg19_gtf_file, config)

# Functions to test
test_that('tb_to_gr_by_position', {
    tb_pos <- tibble(
        pos_col = c('chr1:100-20000', 'chrX:30000-40000', '4:300-500'),
        metacol1 = c(1, 2, 3),
        metacol2 = list(c('A', 'B'), c(1, 2), c(5))
    )
    gr_pos <- tibble(
        seqnames = c('1', 'X', '4'),
        start = c(100, 30000, 300),
        end = c(20000, 40000, 500),
        pos_col = c('chr1:100-20000', 'chrX:30000-40000', '4:300-500'),
        metacol1 = c(1, 2, 3),
        metacol2 = list(c('A', 'B'), c(1, 2), c(5))
    ) %>% as_granges()

    expect_equal(tb_to_gr_by_position(tb_pos, 'pos_col'), gr_pos)

    tb_pos <- bind_rows(
        tb_pos,
        tibble(pos_col = 'XXX:1234')
    )
    expect_error(tb_to_gr_by_position(tb_pos, 'pos_col'))

    #test empty input; expect_qual is tricky due to info in seqnames
    expect_length(
        tb_to_gr_by_position(tb_pos %>% filter(metacol1 == 5), 'pos_col'),
        0
    )
})

#Note: somehow this didn't/doesn't catch all possible issues
# function failed before, due to not correctly checking for all matched gband names
test_that('tb_to_gr_by_gband', {
    
    tb_band <- tibble(
        band_col = c('1p35', '1p11', '1q21.2'),
        metacol1 = c(1, 2, 3),
        metacol2 = list(c('A', 'B'), c(1, 2), c(5))
    )
    gr_band <- tibble(
        seqnames = c('chr1', 'chr1', 'chr1'),
        start = c(28000000, 120600000, 147000000),
        end = c(34600000, 125000000, 150300000),
        band_col = c('1p35', '1p11', '1q21.2'),
        metacol1 = c(1, 2, 3),
        metacol2 = list(c('A', 'B'), c(1, 2), c(5))
    ) %>% as_granges()

    expect_equal(tb_to_gr_by_gband(tb_band,  gr_info, 'band_col'), gr_band)

    # Test error on wrong formatted gband
    tb_band2 <- bind_rows(
        tb_band,
        tibble(band_col = 'XXX')
    )
    expect_error(tb_to_gr_by_gband(tb_band2,  gr_info, 'band_col'))

    # Test error on non-annotated gband
    tb_band2 <- bind_rows(
        tb_band,
        tibble(band_col = '1q44')
    )
    expect_error(tb_to_gr_by_gband(tb_band2,  gr_info, 'band_col'))

    #test empty input; expect_qual is tricky due to info in seqnames
    expect_length(
        tb_to_gr_by_gband(tb_band %>% filter(metacol1 == 5), gr_info, 'band_col'),
        0
    )

})

test_that('parse_hotspot_table', {
    tb <- load_hotspot_table(config)
    target_chrom_style <- 'UCSC'

    expected_gr <- minimal_probes %>%
        mutate(
            seqnames = 'chr1',
            start = c(11873, 28050000, 40000, 0, 30200000, 142600000),
            end = c(14409, 28070000, 50000, 7200000, 32400000, 155000000),
            strand = c('+', '+', '*', '*', '*', '*'),
        ) %>%
        as_granges()

    expect_equal(parse_hotspot_table(tb, gr_genes, gr_info, target_chrom_style), expected_gr)
    
    tb <- bind_rows(tb, tibble(mapping = 'unclear'))
    expect_error(parse_hotspot_table(tb, gr_genes, gr_info, target_chrom_style))
    
    #test empty input
    expect_length(
        parse_hotspot_table(tb %>% filter(description_doi == '123'), gr_genes, gr_info, target_chrom_style),
        0
    )
})

# This requires the default cache files for hg19, do *not* run this in CI
# Also parsing the whole hg19 gtf takes a bit, so allow manual skipping as well?
test_that('parse inbuilt tables', {
    skip_on_ci()
    skip_on_covr()
    config <- list(
        'genome_version' = 'hg19',
        'snakedir' = test_path('../../stemcnv_check/'),
        'global_settings' = list(
            'hg19_gtf_file' = '~/.cache/stemcnv-check/static-data/gencode.hg19.v45.gtf.gz',
            'hg19_genomeInfo_file' = '~/.cache/stemcnv-check/static-data/UCSC_hg19_chromosome-info.tsv'
        ),
        'settings' = list(
            'CNV_processing' = list(
                'gene_overlap' = list(
                    'exclude_gene_type_regex' = c(),
                    'include_only_these_gene_types' = c('lncRNA', 'miRNA', 'protein_coding'),
                    'stemcell_hotspot_list' = '__inbuilt__/supplemental-files/genelist-stemcell-hotspots.tsv',
                    'cancer_gene_list' = '__inbuilt__/supplemental-files/genelist-cancer-drivers.tsv'
                )
            ),
            'SNV_analysis' = list(
                'snv_hotspot_table' = '__inbuilt__/supplemental-files/SNV-stemcell-hotspots.tsv'
            )
        )
    )
    gtf_file <- config$global_settings$hg19_gtf_file
    ginfo_file <- config$global_settings$hg19_genomeInfo_file

    stemcell_hotspot_tb <- load_hotspot_table(config)
    cancer_gene_tb <- load_hotspot_table(config, 'cancer_gene')
    config$settings$CNV_processing$gene_overlap$cancer_gene_list <- '__inbuilt__/supplemental-files/genelist-cancer-hotspots.tsv'
    cancer_gene_tb2 <- load_hotspot_table(config, 'cancer_gene')
    gr_info <- load_genomeInfo(ginfo_file, config)
    gr_genes <- load_gtf_data(gtf_file, config)
    score_settings <- list(
        'pHaplo_threshold' = 0.86,
        'pTriplo_threshold' = 0.94,
        'dosage_sensitive_gene' =  5
    )

    expect_no_error(parse_hotspot_table(stemcell_hotspot_tb, gr_genes, gr_info))
    expect_no_error(parse_hotspot_table(cancer_gene_tb, gr_genes, gr_info))
    expect_no_error(parse_hotspot_table(cancer_gene_tb2, gr_genes, gr_info))
    expect_no_error(get_dosage_sensivity_tb(score_settings) %>% parse_hotspot_table(gr_genes, gr_info))
    # This one can not be parsed like the others (it's never used as a gr either)
    expect_no_error(load_hotspot_table(config, 'snv_hotspot'))
})

test_that('load_hotspot_table', {
    load_hotspot_table(config, 'stemcell_hotspot') %>%
        # remove 'spec_tbl_df' class from readr
        .[] %>% expect_equal(minimal_probes)
    load_hotspot_table(config, 'cancer_gene') %>%
        .[] %>% expect_equal(minimal_probes)
})

# get_roi_gr <- function(sample_id, sampletable, config, gr_genes, gr_info, target_chrom_style)
test_that('get_roi_gr', {
    
    sample_id <- 'testsample'
    sampletable <- tibble(
        Sample_ID = 'testsample',
        Regions_of_Interest = 'DDX11L1;important|1q21.1;edit site|1:1000-2000;mygene|dummyC'
    )
    target_chrom_style <- get_target_chrom_style(config, gr_genes)
    
    expected <- GRanges(
        seqnames = c('chr1', 'chr1', 'chr1', 'chr1'),
        ranges = IRanges(
            start = c(11873, 142600000, 1000, 28050000),
            end = c(14409, 147000000, 2000, 28070000)
        ),
        strand = c('+', '*', '*', '+'), 
        list_name = 'ROI',
        hotspot = c('DDX11L1', '1q21.1', '1:1000-2000', 'dummyC'),
        mapping = c('gene_name', 'gband', 'position', 'gene_name'),
        call_type = 'any',
        check_score = 50,
        description = c('ROI_1: DDX11L1', 'ROI_2: 1q21.1; important', 'ROI_3: 1:1000-2000; edit site', 'ROI_4: dummyC; mygene'),
        description_doi = NA_character_,
        description_htmllinks = c('ROI_1: DDX11L1', 'ROI_2: 1q21.1; important', 'ROI_3: 1:1000-2000; edit site', 'ROI_4: dummyC; mygene')
    )
    
    get_roi_gr(sample_id, sampletable, config, gr_genes, gr_info, target_chrom_style) %>%
        expect_equal(expected)
    
    sampletable$Regions_of_Interest <- NA_character_
    empty <- GRanges(
        list_name = character(),
        hotspot = character(),
        mapping = character(),
        call_type = character(),
        check_score = numeric(),
        description = character(),
        description_doi = character(),
        description_htmllinks = character()
    )
    get_roi_gr(sample_id, sampletable, config, gr_genes, gr_info, target_chrom_style) %>%
        expect_equal(empty)
})