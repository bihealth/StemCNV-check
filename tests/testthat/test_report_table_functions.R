library(tidyverse)
library(plyranges)
library(DT)
library(knitr)
library(testthat)

source(test_path("../../stemcnv_check/scripts/R/helper_functions.R"))
source(test_path("../../stemcnv_check/scripts/R/report_table_functions.R"))

# Functions to test:
# - [ ] vector_to_js
# - [ ] format_column_names
# - [ ] simple_table_output
# - [ ] summary_table
# - [x] format_hotspots_to_badge
# - [ ] CNV_table_output
# - [ ] gene_table_output
# - [x] hotspot_table_output
# - [ ] SNV_table_output

config <- list(
    'snakedir' = '',
    'settings' = list(
        'CNV_processing' = list(
            'gene_overlap' = list(
                'stemcell_hotspot_list' = test_path('../data/minimal-hotspots.tsv'),
                'cancer_gene_list' = test_path('../data/minimal-hotspots.tsv')
            )
        )
    )
)

# Test `format_hotspots_to_badge` function
# format_hotspots_to_badge <- function(hotspot_vec, CNVtype_vec, gene_details, color = 'red')
test_that("format_hotspots_to_badge", {
    testthat::local_edition(3)
    hotspot_vec <- c("", "1q21", "1q21", "dummyC", "1p36|DDX11L1", "1p36|DDX11L1")
    CNVtype_vec <- c("gain", "gain", "LOH", "gain", "loss", "gain")
    # 1 - empty
    # 2 - gband hit, matching CNV (gain)
    # 3 - gband hit, not matching CNV (Note: this case should never occur; warning needed?)
    # 4 - gene hit
    # 5 - gene hit & gband hit matching CNV (loss)
    # 6 - gene hit & gband hit not matching CNV (Note: theoretically possible)

    gene_details <- load_hotspot_table(config, 'stemcell_hotspot')
    
    expected <- c(
        '-', 
        '<span class="badge badge-red" title="test-list&#013;Check_Score contribution: 10&#013;Sources: dummy{1},dummy{2}">1q21</span>', 
        '1q21', 
        '<span class="badge badge-red" title="test-list&#013;Check_Score contribution: 15&#013;Something: Dummy{1}">dummyC</span>', 
        paste0(
            '<span class="badge badge-red" title="test-list&#013;Check_Score contribution: 10&#013;',
            'Sources: dummy{1}&#013;Something: else{2}">1p36</span>',
            '<span class="badge badge-red" title="test-list&#013;Check_Score contribution: 30&#013;',
            'Sources: dummy">DDX11L1</span>'
        ),
        '1p36<span class="badge badge-red" title="test-list&#013;Check_Score contribution: 30&#013;Sources: dummy">DDX11L1</span>'
    )
    expect_equal(
        format_hotspots_to_badge(hotspot_vec, CNVtype_vec, gene_details, 'red'),
        expected
    )
    
    #test with include_hover = FALSE & shorthand = orange
    expected <- c(
        '-', 
        '<span class="badge badge-orange">1q21</span>', 
        '1q21', 
        '<span class="badge badge-orange">dummyC</span>', 
        '<span class="badge badge-orange">1p36</span><span class="badge badge-orange">DDX11L1</span>', 
        '1p36<span class="badge badge-orange">DDX11L1</span>'
    )
    expect_equal(
        format_hotspots_to_badge(hotspot_vec, CNVtype_vec, gene_details, 'orange', FALSE),
        expected
    )
})

# hotspot_table_output(hotspots, cnv_type, plotsection, stemcell_hotspot_tb, cancer_gene_tb, report_config, out_format) %>%
test_that("hotspot_table_output", {
    hotspots <- c('DDX11L1', '1p36')
    cnv_type <- 'loss'
    stemcell_hotspot_tb <- load_hotspot_table(config, 'stemcell_hotspot') 
    dosage_sensitive_gene_tb <- tibble(
        list_name = NA_character_,
        hotspot = NA_character_,
        call_type = NA_character_,
        description = NA_character_,
        check_score = NA_real_,
        mapping = NA_character_,
        description_doi = NA_character_
    )
    cancer_gene_tb <- dosage_sensitive_gene_tb
    # these aren't used so far
    plotsection <- 'test'
    report_config <- list()
    
    expected_tb <- stemcell_hotspot_tb %>%
        filter(hotspot %in% hotspots) 
    
    expected <- expected_tb %>%
        select(hotspot, call_type, list_name, description_htmllinks, check_score, mapping, description_doi) %>%
        dplyr::rename(description = description_htmllinks) %>%
        rename_with(format_column_names) %>%
        mutate(
            Hotspot = paste0('<span class="badge badge-red">', Hotspot, '</span>'),
            Description = str_replace_all(Description, '&#013;', '<br/>')
        ) %>%
        datatable(
            options = list(
                dom = 't', 
                #extensions = c('Buttons'),
                buttons = c('colvis', 'copy', 'print'),
                pageLength = 2,
                columnDefs = list(list(targets = 5:6, visible = FALSE))
            ),
            rownames = FALSE,
            escape = FALSE
        )
    
    hotspot_table_output(
        hotspots, cnv_type, plotsection, 
        stemcell_hotspot_tb, dosage_sensitive_gene_tb, cancer_gene_tb,
        report_config, 'html'
    ) %>%
        expect_equal(expected)
    
    # test non-html output
    expected <- expected_tb %>%
        select(hotspot, call_type, list_name, description, check_score, description_doi) %>%
        dplyr::rename(dois = description_doi) %>%
        rename_with(format_column_names) %>% 
        kable()
    hotspot_table_output(
        hotspots, cnv_type, plotsection, 
        stemcell_hotspot_tb, dosage_sensitive_gene_tb, cancer_gene_tb,
        report_config, 'not-html'
    ) %>%
        expect_equal(expected)
    
    # test with cancer_gene table
    cancer_gene_tb <- stemcell_hotspot_tb %>%
        filter(hotspot == '1p36') %>%
        mutate(list_name = 'cancer_gene')
    stemcell_hotspot_tb.no_ov <- stemcell_hotspot_tb %>%
        filter(hotspot != '1p36')
    expected <- expected_tb %>% 
        mutate(list_name = ifelse(hotspot == '1p36', 'cancer_gene', list_name)) %>%
        select(hotspot, call_type, list_name, description, check_score, description_doi) %>%
        dplyr::rename(dois = description_doi) %>%
        rename_with(format_column_names) %>% 
        kable()
    hotspot_table_output(
        hotspots, cnv_type, plotsection,
        stemcell_hotspot_tb.no_ov, dosage_sensitive_gene_tb, cancer_gene_tb,
        report_config, 'not-html'
    ) %>%
        expect_equal(expected)
    # test with same hotspot in both stemcell_hotspot and cancer_gene table
    expected <- expected_tb %>% 
        mutate(list_name = ifelse(hotspot == '1p36', 'test-list|cancer_gene', list_name)) %>%
        separate_rows(list_name, sep = '\\|') %>%
        select(hotspot, call_type, list_name, description_htmllinks, check_score, mapping, description_doi) %>%
        dplyr::rename(description = description_htmllinks) %>%
        rename_with(format_column_names) %>%
        mutate(
            Hotspot = paste0('<span class="badge badge-', c('red', 'red', 'orange'), '">', Hotspot, '</span>'),
            Description = str_replace_all(Description, '&#013;', '<br/>')
        ) %>%
        datatable(
            options = list(
                dom = 't', 
                # extensions = c('Buttons'),
                buttons = c('colvis', 'copy', 'print'),
                pageLength = 3,
                columnDefs = list(list(targets = 5:6, visible = FALSE))
            ),
            rownames = FALSE,
            escape = FALSE
        )
    hotspot_table_output(
        hotspots, cnv_type, plotsection, 
        stemcell_hotspot_tb, dosage_sensitive_gene_tb, cancer_gene_tb,
        report_config, 'html'
    ) %>%
        expect_equal(expected)
    # test with only partially matching cnv_type
    cnv_type <-  'LOH'
    expected <- expected_tb %>%
        filter(hotspot == 'DDX11L1') %>%
        select(hotspot, call_type, list_name, description, check_score, description_doi) %>%
        dplyr::rename(dois = description_doi) %>%
        rename_with(format_column_names) %>% 
        kable()
    hotspot_table_output(
        hotspots, cnv_type, plotsection,
        stemcell_hotspot_tb, dosage_sensitive_gene_tb, cancer_gene_tb,
        report_config, 'not-html'
    ) %>%
        expect_equal(expected)
    
})