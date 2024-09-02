library(testthat)

library(tidyverse)
library(plyranges)

source(test_path('../../stemcnv_check/scripts/R/helper_functions.R'))

# Functions to test:
# - [x] read_sampletable
# - [x] get_sample_info
# - [x] fix_CHROM_format
# - [x] get_sex_chroms
# - [x] get_target_chrom_style
# - [ ] load_preprocessed_cnvs  #! consider replacing this with read_CNV_vcf (+ annotation support)
# - [ ] get_static_path         #! hard to test; should only do things if inside the Rmd environment; deprecate?
# - [x] load_gtf_data
# - [x] load_genomeInfo
# - [ ] ensure_list_cols
# - [ ] get_SNP_clustering_IDs  #! this should be replaced by direct snakemake handover 

# Note: primary checks are done in python, R just needs to work
test_that("read_sampletable", {
    expected_tb <- tibble(
        Sample_Name = c("Cellline-A MasterBank", "Cellline-A WorkingBank", "Cellline-B MasterBank", "Cellline-B-1 clone1"),
        Chip_Name = c("123456789000", "123456789001", "123456789000", "123456789001"),
        Chip_Pos = c("R01C01", "R01C01", "R01C02", "R01C02"),
        Sample_ID = c("Cellline-A-MB", "Cellline-A-WB", "Cellline-B-MB", "Cellline-B-1-cl1"),
        Sex = c("Female", "Female", "Male", "Male"),
        Reference_Sample = c(NA, "Cellline-A-MB", NA, "Cellline-B-MB"),
        Sample_Group = c("ExampleCellines", "ExampleCellines", "ExampleCellines", "ExampleCellines"),
        Regions_of_Interest = c(
                "Example1|chr1:100000-200000", 
                "Example1|chr1:100000-200000;chr11:60000-70000", 
                "Example1|chr1:100000-200000;chr11:60000-70000",
                NA
        )
    )
    
    test_path('../../stemcnv_check/control_files/sample_table_example.tsv') %>%
        read_sampletable() %>%
        # this removes the "spec_tbl_df" class that readr adds
        # see here: https://www.tidyverse.org/blog/2018/12/readr-1-3-1/
        .[] %>%
        expect_equal(expected_tb)
})

test_that("get_sample_info", {
    sampletable <- test_path('../../stemcnv_check/control_files/sample_table_example.tsv') %>%
        read_sampletable()
    
    expect_equal(get_sample_info('Cellline-A-MB', 'ref_id', sampletable), NA_character_)
    expect_equal(get_sample_info('Cellline-A-MB', 'sex', sampletable), 'f')
    
    expect_equal(get_sample_info('Cellline-B-1-cl1', 'ref_id', sampletable), 'Cellline-B-MB')
    expect_equal(get_sample_info('Cellline-B-1-cl1', 'sex.ref', sampletable), 'm')
    expect_equal(get_sample_info('Cellline-B-MB', 'sex', sampletable), 'm')
    
    expect_error(get_sample_info('Cellline-B-1-cl1', 'unsupported', sampletable))
    sampletable$Sex <- c('f', 'f', 'f', 'm')
    expect_error(get_sample_info('Cellline-B-1-cl1', 'sex.ref', sampletable))
})

input_tb <- tibble(
    seqnames = c('2', 'X', '1'),
    start = c(1e4, 3e6, 53e5),
    end = c(3e4, 5e6, 77e5)
)

# there is also a test for this is test_preprocess_CNV_functions.R
test_that("fix_CHROM_format", {
    expected_ncbi <- input_tb %>%
        mutate(seqnames = factor(c('2', 'X', '1'), levels = c('1', '2', 'X'))) %>%
        as_granges()
    expected_ucsc <- input_tb %>%
        mutate(seqnames = factor(c('chr2', 'chrX', 'chr1'), levels = c('chr1', 'chr2', 'chrX'))) %>%
        as_granges()
    
    expect_equal(as_granges(input_tb) %>% fix_CHROM_format('NCBI'), expected_ncbi)
    expect_equal(as_granges(input_tb) %>% fix_CHROM_format('UCSC'), expected_ucsc)
})


test_that("get_sex_chroms", {
    expect_equal(get_sex_chroms(input_tb), c('X', 'Y'))
    expect_equal(get_sex_chroms(as_granges(input_tb)), c('X', 'Y'))
    as_granges(input_tb) %>%
        fix_CHROM_format('UCSC') %>%
        get_sex_chroms() %>%
        expect_equal(c('chrX', 'chrY'))
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