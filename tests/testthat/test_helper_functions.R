library(testthat)

library(tidyverse)
library(plyranges)

source(test_path('../../src/stemcnv_check/scripts/R/helper_functions.R'))

# Functions to test:
# - [x] read_sampletable
# - [x] get_sample_info
# - [x] fix_CHROM_format
# - [x] get_sex_chroms
# - [x] get_target_chrom_style
# - [ ] fix_rel_filepath         #! hard to test; should only do things if inside the Rmd environment; deprecate?
# - [x] load_gtf_data
# - [x] load_genomeInfo
# - [x] load_hotspot_table
# - [x] split_merged_CNV_callers

# Note: primary consistency checks are done in python, R just needs to work
test_that("read_sampletable", {
    expected_tb <- tibble(
        Sample_ID = c("Cellline-A-MB", "Cellline-A-WB", "Cellline-B-MB", "Cellline-B-1-cl1"),
        Chip_Name = c("123456789000", "123456789001", "123456789000", "123456789001"),
        Chip_Pos = c("R01C01", "R01C01", "R01C02", "R01C02"),
        Array_Name = rep('ExampleArray', 4),
        Sex = c("Female", "Female", "Male", "Male"),
        Reference_Sample = c(NA, "Cellline-A-MB", NA, "Cellline-B-MB"),
        Regions_of_Interest = c(
                "Example1|chr1:100000-200000", 
                "Example1|chr1:100000-200000;chr11:60000-70000", 
                "Example1|chr1:100000-200000;chr11:60000-70000",
                NA
        ),
        Sample_Group = c("ExampleCellines", "ExampleCellines", "ExampleCellines", "ExampleCellines"),
    )
    
    test_path('../../src/stemcnv_check/control_files/sample_table_example.tsv') %>%
        read_sampletable() %>%
        # this removes the "spec_tbl_df" class that readr adds
        # see here: https://www.tidyverse.org/blog/2018/12/readr-1-3-1/
        .[] %>%
        expect_equal(expected_tb)
    
    # test that regex name change works, use tmp file instead of fake filesystem
    tmp_file <- tempfile(fileext='.tsv')
    expected_tb %>%
        rename_with(~paste(., 'dummy')) %>%
        write_tsv(tmp_file)
    read_sampletable(tmp_file, col_remove_regex = ' .*') %>% expect_equal(expected_tb)
})

test_that("get_sample_info", {
    sampletable <- test_path('../../src/stemcnv_check/control_files/sample_table_example.tsv') %>%
        read_sampletable()
    config <- list()
    
    expect_equal(get_sample_info('Cellline-A-MB', 'ref_id', config, sampletable), NA_character_)
    expect_equal(get_sample_info('Cellline-A-MB', 'sex', config, sampletable), 'f')
    
    expect_equal(get_sample_info('Cellline-B-1-cl1', 'ref_id', config, sampletable), 'Cellline-B-MB')
    expect_equal(get_sample_info('Cellline-B-MB', 'sex', config, sampletable), 'm')
    
    expect_error(get_sample_info('Cellline-B-1-cl1', 'unsupported', config, sampletable))
})

input_tb <- tibble(
    seqnames = c('2', 'X', '1'),
    start = c(1e4, 3e6, 53e5),
    end = c(3e4, 5e6, 77e5)
)

# there is also a test for this in test_CNV_preprocess_functions.R
test_that("fix_CHROM_format", {
    expected_ncbi <- input_tb %>%
        mutate(seqnames = factor(c('2', 'X', '1'), levels = c('1', '2', 'X')))
    expected_ucsc <- input_tb %>%
        mutate(seqnames = factor(c('chr2', 'chrX', 'chr1'), levels = c('chr1', 'chr2', 'chrX')))
    
    # Test on granges input
    expect_equal(as_granges(input_tb) %>% fix_CHROM_format('NCBI'), as_granges(expected_ncbi))
    expect_equal(as_granges(input_tb) %>% fix_CHROM_format('UCSC'), as_granges(expected_ucsc))
    
    # Test on tibble input
    expect_equal(input_tb %>% fix_CHROM_format('NCBI'), expected_ncbi)
    expect_equal(input_tb %>% fix_CHROM_format('UCSC'), expected_ucsc)
    
    # Expect error on invalid input
    expect_error(input_tb %>% fix_CHROM_format('invalid'))
})


test_that("get_sex_chroms", {
    expect_equal(get_sex_chroms('NCBI'), c('X', 'Y'))
    expect_equal(get_sex_chroms('UCSC'), c('chrX', 'chrY'))
    expect_error(get_sex_chroms('invalid'))
})


test_that('get_target_chrom_style', {
    config <- list(
        settings = list(
            vcf_output = list(
                chrom_style = 'keep-original'
            )
        )
    )
    gr <- as_granges(input_tb)
    
    expect_equal(get_target_chrom_style(config, gr), 'NCBI')
    config$settings$vcf_output$chrom_style <- 'UCSC'
    expect_equal(get_target_chrom_style(config, gr), 'UCSC')
})

test_that("load_gtf_data", {
    # Construct minimal config
    config <- list(
        'genome_version' = 'hg19',
        'global_settings' = list(
            'hg19_gtf_file' = test_path('../data/hg_minimal.gtf')
        ),
        'settings' = list(
            'CNV_processing' = list(
                'gene_overlap' = list(
                    'exclude_gene_type_regex' = c(),
                    'include_only_these_gene_types' = c('lncRNA', 'miRNA', 'protein_coding'),
                    'whitelist_hotspot_genes' = FALSE
                )
            ),
            'vcf_output' = list(
                'chrom_style' = 'UCSC'
            )
        )
    )
    gtf_file <- test_path('../data/hg_minimal.gtf')
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

    # Test default gene selection ('lncRNA', 'miRNA', 'protein_coding')
    expect_equal(
        load_gtf_data(gtf_file, config),
        as_granges(expected_tb %>% filter(gene_type %!in% c("unprocessed_pseudogene", 'something_else')))
    )
    # Test additional gene type included, and another one excluded
    config$settings$CNV_processing$gene_overlap$include_only_these_gene_types <- c(
        'lncRNA', 'miRNA', 'protein_coding', 'something_else'
    )
    config$settings$CNV_processing$gene_overlap$exclude_gene_type_regex <- c('lncRNA')
    expect_equal(
        load_gtf_data(gtf_file, config),
        as_granges(expected_tb %>% filter(gene_type %!in% c('lncRNA', 'unprocessed_pseudogene')))
    )
    # Test selection only by exclusion, no specific included types
    config$settings$CNV_processing$gene_overlap$include_only_these_gene_types <- c()
    config$settings$CNV_processing$gene_overlap$exclude_gene_type_regex <- c(
        'lncRNA', 'something_else', 'unprocessed_pseudogene'
    )
    expect_equal(load_gtf_data(gtf_file, config), as_granges(expected_tb %>% filter(gene_type == 'protein_coding')))
    # Test no gene type filtering
    config$settings$CNV_processing$gene_overlap$include_only_these_gene_types <- c()
    config$settings$CNV_processing$gene_overlap$exclude_gene_type_regex <- c()
    expect_equal(load_gtf_data(gtf_file, config), as_granges(expected_tb))
    # Test inclusion of specific gene_names
    config$settings$CNV_processing$gene_overlap$whitelist_hotspot_genes <- TRUE
    config$settings$CNV_processing$gene_overlap$include_only_these_gene_types <- c('lncRNA', 'miRNA', 'protein_coding')
    expect_equal(
        load_gtf_data(gtf_file, config, include_hotspot_genes = c('dummyB')),
        as_granges(expected_tb %>% filter(gene_type %!in% c("unprocessed_pseudogene")))
    )    
})

test_that('load_genomeInfo', {
    config <- list(
        'genome_version' = 'hg19',
        'global_settings' = list(
            'hg19_genomeInfo_file' = test_path('../data/gr_info_minimal.tsv')
        )
    )
    ginfo_file <- test_path('../data/gr_info_minimal.tsv')
    
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

    expect_equal(load_genomeInfo(ginfo_file, config), expected_gr)
})

test_that('split_merged_CNV_callers', {
    # Example calls defined in helper_combined_cnv_data.R
    # source(test_path('helper_combined_cnv_data.R'))
    
    #Note: only some value are restored for initial_call_details string:
    # <caller>_<start>-<end>_CN<copy_number>_cov<coverage_merge_call>_<filter>
    # other values (i.e. n_probes, probe_density_Mb) are kept from the recalculation for the merged call
    expected <- bind_rows(toolA, toolB) %>%
        as_granges() %>%
        # expect same order as input
        select(colnames(combined_tools@elementMetadata)) %>%
        # will not be the same
        select(sample_id, CNV_type, CNV_caller, CN, FILTER) %>%
        arrange(CNV_caller, start)
    
    
    split_merged_CNV_callers(combined_tools, defined_labels) %>%
        select(sample_id, CNV_type, CNV_caller, CN, FILTER) %>%
        expect_equal(expected)
})