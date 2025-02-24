library(tidyverse)
library(plyranges)
library(DT)
library(knitr)
library(testthat)

source(test_path("../../stemcnv_check/scripts/R/helper_functions.R"))
source(test_path("../../stemcnv_check/scripts/R/report_table_functions.R"))

defined_labels <- yaml.load_file(test_path('../../stemcnv_check/control_files/label_name_definitions.yaml'))

# Functions to test:
# - [ ] vector_to_js
# - [ ] format_column_names
# - [ ] simple_table_output
# - [x] summary_table
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
    ),
    'evaluation_settings' = list(
        'summary_stat_warning_levels' = list(
            'call_count_excl_labels' = list('Excluded'),
            'use_last_level' = list('call_rate', 'computed_gender', 'SNP_pairwise_distance_to_reference')
        )
    )
)

# Test `format_hotspots_to_badge` function
# format_hotspots_to_badge <- function(hotspot_vec, CNVtype_vec, gene_details, color = 'red')
test_that("format_hotspots_to_badge", {
    testthat::local_edition(3)
    hotspot_vec <- c("", "1q21", "1q21", "dummyC", "1p36|DDX11L1", "1p36|DDX11L1", 'A|B|C')
    CNVtype_vec <- c("gain", "gain", "LOH", "gain", "loss", "gain", 'gain')
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
        '1p36<span class="badge badge-red" title="test-list&#013;Check_Score contribution: 30&#013;Sources: dummy">DDX11L1</span>',
        'A; B; C'
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
        '1p36<span class="badge badge-orange">DDX11L1</span>',
        'A; B; C'
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

#summary_table(summary_stat_table, sample_headers, config, defined_labels)
test_that("summary_table", {
    
    summary_stat_table <- tibble(
        Description = c(
            'sample_id', 'call_rate', 'computed_gender', 'SNP_pairwise_distance_to_reference',
            'SNPs_post_filter', defined_labels$sample_qc_measures[3:10]
        ),
        sample_value = c(
            'SampleID', '0.991', 'M', '60% (123 SNPs)', '123456', '-0.5', '50', '30', '15', '7', '0', '0', '0'
        ),
        sample_eval = c(
            'SampleID', 'OK', 'OK', NA, 'high concern', 'OK', 'warning', 'unusual', 'warning', 'unusual', 'OK', 'OK', 'OK'
        ),
        reference_value = c(
            'SampleID2', '0.995', 'M', '61% (135 SNPs)', NA, '0.3', '40', '20', NA, NA, NA, NA, NA
        ),
        reference_eval = c(
            'SampleID2', 'OK', 'OK', NA, NA, 'OK', 'unusual', 'OK', NA, NA, NA, NA, NA
        ),
    )
    
    sample_headers <- set_names(
        c('SampleID', 'Reference (SampleID2)'),
        c('SampleID', 'SampleID2')
    )
    
    expected_tb <- summary_stat_table %>%
        filter(Description != 'sample_id') %>%
        select(-contains('eval')) %>%
        mutate(Description = format_column_names(Description)) %>%
        set_names(c(' ', sample_headers))
    
    green <- 'rgb(146,208,80)'
    expected_colors <- tibble(
        SampleID = c(green, green, 'white', 'red', green, 'orange', 'yellow', 'orange', 'yellow', green, green, green),
        SampleID2 = c(green, green, 'white', 'white', green, 'yellow', green, 'white', 'white', 'white', 'white', 'white'),
    ) %>% set_names(sample_headers)
    
    ignored_calls <- paste0(
        '\\nCalls with one of these Labels are not counted: ',
        config$evaluation_settings$summary_stat_warning_levels$call_count_excl_labels %>% paste(collapse = '|')
    )
    summary_row_help <- c(
        "Call Rate" = paste0(
            'The Illumina call rate corresponds to the percentage of probes for which a clear genotype could be ',
            'determined.\\nLow values are strong indicator of sample quality issues that also impact make any ',
            'further analysis including CNV calling.'
        ),
        "Computed Gender" = paste0(
            'The Illumina computed gender (sex) based on X and Y chromosome probes.\\n',
            'Mismatches with annotated sex can indicate annotation mistakes, sample swaps, or severe quality issues.'
        ),
        "SNPs Post Filter" = 'The percentage of SNP probes retained after the employed StemCNV-check filter strategy.',
        "SNP Pairwise Distance To Reference" = paste0(
            'The number of probes with a different genotype than the reference sample.\\n',
            'Calculated as pairwise difference, which may include more probes than used for sample clustering.\\n',
            'Increased values can indicate a sample swap or a considerable number of genomic mutations between ',
            'sample and references.'
        ),
        "Loss Gain Log2ratio" = paste0(
            'The log2 transformed ratio of loss and gain CNV calls.\\n',
            'Deviation from equal balance (0) can indicate potential quality issues or problems with CNV calling.',
            ignored_calls
        ),
        "Total Calls CNV" = paste0('The total number of CNV (gain/loss) calls.', ignored_calls),
        "Total Calls LOH" = paste0('The total number of LOH calls.', ignored_calls),
        "Reportable Calls CNV" = 'The number of CNV calls designated as "reportable".',
        "Reportable Calls LOH" = 'The number of LOH calls designated as "reportable".',
        "Critical Calls CNV" = 'The number of CNV calls designated as "critical".',
        "Critical Calls LOH" = 'The number of LOH calls designated as "critical".',
        "Critical SNVs" = 'The number of detected SNVs designated as "critical".'
    )
    
    expected <- expected_tb %>%
        datatable(
            options = list(
                dom = 't',
                pageLength = nrow(expected_tb),
                rowCallback = JS(
                    "function(row, data, displayNum, displayIndex, dataIndex) {",
                    "let help_text = ", vector_to_js(summary_row_help), ";",
                    #"console.log('hover-test: ' + help_text[data[0]]);",
                    "$('td', row).attr('title', help_text[data[0]]);",
                    "}"
                )
            ),
            rownames = FALSE
        ) %>% 
        # Color coding of values
        formatStyle(
            2, 
            backgroundColor = styleRow(1:nrow(expected_tb), unlist(expected_colors[, sample_headers[[1]]])),
            textAlign = 'center'
        ) %>%
        formatStyle(
            3,
            backgroundColor = styleRow(1:nrow(expected_tb), unlist(expected_colors[, sample_headers[[2]]])), 
            textAlign = 'center'
        ) %>%
        formatStyle(
            1, 
            fontWeight = styleEqual(
                unlist(config$evaluation_settings$summary_stat_warning_levels$use_last_level) %>% format_column_names(), 
                'bold'
            )
        )
    
    expect_equal(
        summary_table(summary_stat_table, sample_headers, config, defined_labels),
        expected
    )
    
    # also test without reference
    summary_stat_table <- summary_stat_table %>%
        select(-contains('reference'))
    sample_headers <- sample_headers[1]
    
    expected <- expected_tb %>%
        select(-3) %>%
        datatable(
            options = list(
                dom = 't',
                pageLength = nrow(expected_tb),
                rowCallback = JS(
                    "function(row, data, displayNum, displayIndex, dataIndex) {",
                    "let help_text = ", vector_to_js(summary_row_help), ";",
                    #"console.log('hover-test: ' + help_text[data[0]]);",
                    "$('td', row).attr('title', help_text[data[0]]);",
                    "}"
                )
            ),
            rownames = FALSE
        ) %>% 
        # Color coding of values
        formatStyle(
            2, 
            backgroundColor = styleRow(1:nrow(expected_tb), unlist(expected_colors[, sample_headers[[1]]])),
            textAlign = 'center'
        ) %>%
        formatStyle(
            3, 
            backgroundColor = styleRow(1:nrow(expected_tb), unlist(expected_colors[, sample_headers[[1]]])),
            textAlign = 'center'
        ) %>%
        # Make all last level (= potentially red) rows have bold text
        formatStyle(
            1, 
            fontWeight = styleEqual(
                unlist(config$evaluation_settings$summary_stat_warning_levels$use_last_level) %>% format_column_names(), 
                'bold'
            )
        )
    
    expect_equal(
        summary_table(summary_stat_table, sample_headers, config, defined_labels),
        expected
    )    
})
