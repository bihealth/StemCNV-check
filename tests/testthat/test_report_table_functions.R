library(tidyverse)
library(plyranges)
library(DT)
library(knitr)
library(testthat)

source(test_path("../../src/stemcnv_check/scripts/R/helper_functions.R"))
source(test_path("../../src/stemcnv_check/scripts/R/hotspot_functions.R"))
source(test_path("../../src/stemcnv_check/scripts/R/report_table_functions.R"))

# Functions to test:
# - [ ] vector_to_js
# - [ ] format_column_names
# - [ ] simple_table_output
# - [x] summary_table
# - [x] format_hotspots_to_badge
# - [ ] CNV_table_output
# - [ ] gene_table_output
# - [x] hotspot_table_output
# - [x] SNV_table_output

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
                'stemcell_hotspot_list' = test_path('../data/minimal-hotspots.tsv'),
                'cancer_gene_list' = test_path('../data/minimal-hotspots.tsv'),
                'whitelist_hotspot_genes' = FALSE,
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
gr_info <- load_genomeInfo(config$global_settings$hg19_genomeInfo_file, config)
gr_genes <- load_gtf_data(config$global_settings$hg19_gtf_file, config)

# format_hotspots_to_badge <- function(hotspot_vec, CNVtype_vec, colorgene_details, gene_details, include_hover = TRUE)
test_that("format_hotspots_to_badge", {
    testthat::local_edition(3)
    hotspot_vec <- c("", "1q21", "1q21", "dummyC", "1p36|DDX11L1", "1p36|DDX11L1", 'A|B|C')
    CNVtype_vec <- c("gain", "gain", "LOH", "gain", "loss", "gain", 'gain')
    color_vec <- c('red', 'orange', 'red', 'red', 'red', 'orange', 'red')
    # 1 - empty
    # 2 - gband hit, matching CNV (gain)
    # 3 - gband hit, not matching CNV (Note: this case should never occur; warning needed?)
    # 4 - gene hit
    # 5 - gene hit & gband hit matching CNV (loss)
    # 6 - gene hit & gband hit not matching CNV (Note: theoretically possible)

    gene_details <- load_hotspot_table(config, 'stemcell_hotspot')
    
    expected <- c(
        '-', 
        '<span class="badge badge-orange" title="test-list&#013;Check_Score contribution: 10&#013;Sources: dummy{1},dummy{2}">1q21</span>', 
        '1q21', 
        '<span class="badge badge-red" title="test-list&#013;Check_Score contribution: 15&#013;Something: Dummy{1}">dummyC</span>', 
        paste0(
            '<span class="badge badge-red" title="test-list&#013;Check_Score contribution: 10&#013;',
            'Sources: dummy{1}&#013;Something: else{2}">1p36</span>',
            '<span class="badge badge-red" title="test-list&#013;Check_Score contribution: 30&#013;',
            'Sources: dummy">DDX11L1</span>'
        ),
        '1p36<span class="badge badge-orange" title="test-list&#013;Check_Score contribution: 30&#013;Sources: dummy">DDX11L1</span>',
        'A; B; C'
    )
    expect_equal(
        format_hotspots_to_badge(hotspot_vec, CNVtype_vec, color_vec, gene_details),
        expected
    )
    
    #test with include_hover = FALSE & fixed color
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
        format_hotspots_to_badge(hotspot_vec, CNVtype_vec, 'orange', gene_details,FALSE),
        expected
    )
    # test error on wrong color
    expect_error(
        format_hotspots_to_badge(hotspot_vec, CNVtype_vec, 'dummy', gene_details, FALSE)
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
            'sample_id', defined_labels$sample_qc_measures
        ),
        sample_value = c(
            'SampleID', '0.991', 'M', '60% (123 SNPs)', '123456', '-0.5', '50', '30', '15', '7', '3', '0', '0', '0'
        ),
        sample_eval = c(
            'SampleID', 'OK', 'OK', NA, 'high concern', 'OK', 'warning', 'unusual', 'warning', 'unusual', 'unusual', 'OK', 'OK', 'OK'
        ),
        reference_value = c(
            'SampleID2', '0.995', 'M', '61% (135 SNPs)', NA, '0.3', '40', '20', NA, NA, NA, NA, NA, NA
        ),
        reference_eval = c(
            'SampleID2', 'OK', 'OK', NA, NA, 'OK', 'unusual', 'OK', NA, NA, NA, NA, NA, NA
        ),
    )
    
    sample_headers <- set_names(
        c('SampleID', 'Reference (SampleID2)'),
        c('SampleID', 'SampleID2')
    )
    
    expected_tb <- summary_stat_table %>%
        filter(Description != 'sample_id') %>%
        select(-contains('eval')) %>%
        mutate(Description = format_column_names(Description))
    
    green <- 'rgb(146,208,80)'
    expected_colors <- tibble(
        SampleID = c(green, green, 'white', 'red', green, 'orange', 'yellow', 'orange', 'yellow', 'yellow', green, green, green),
        SampleID2 = c(green, green, 'white', 'white', green, 'yellow', green, 'white', 'white', 'white', 'white', 'white', 'white'),
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
        "Reportable SNVs" = 'The number of detected SNVs designated as "reportable".',
        "Critical Calls CNV" = 'The number of CNV calls designated as "critical".',
        "Critical Calls LOH" = 'The number of LOH calls designated as "critical".',
        "Critical SNVs" = 'The number of detected SNVs designated as "critical".'
    )
    
    expected_data <- expected_tb[1:7,] %>%
        set_names(c('Data QC measures', sample_headers)) %>%
        datatable(
            options = list(
                dom = 't',
                pageLength = 7,
                rowCallback = JS(
                    "function(row, data, displayNum, displayIndex, dataIndex) {",
                    "let help_text = ", vector_to_js(summary_row_help[1:7]), ";",
                    "$('td', row).attr('title', help_text[data[0]]);",
                    "}"
                )
            ),
            rownames = FALSE
        ) %>% 
        # Color coding of values
        formatStyle(
            2, 
            backgroundColor = styleRow(1:7, unlist(expected_colors[1:7, sample_headers[[1]]])),
            textAlign = 'center'
        ) %>%
        formatStyle(
            3,
            backgroundColor = styleRow(1:7, unlist(expected_colors[1:7, sample_headers[[2]]])), 
            textAlign = 'center'
        ) %>%
        formatStyle(
            1, 
            fontWeight = styleEqual(
                unlist(config$evaluation_settings$summary_stat_warning_levels$use_last_level) %>% format_column_names(), 
                'bold'
            )
        )
    
    expected_sample <- expected_tb[8:13, 1:2] %>%
        set_names(c('Sample QC measures', sample_headers[1])) %>%
        datatable(
            options = list(
                dom = 't',
                pageLength = 6,
                rowCallback = JS(
                    "function(row, data, displayNum, displayIndex, dataIndex) {",
                    "let help_text = ", vector_to_js(summary_row_help[8:13]), ";",
                    "$('td', row).attr('title', help_text[data[0]]);",
                    "}"
                )
            ),
            rownames = FALSE
        ) %>% 
        # Color coding of values
        formatStyle(
            2, 
            backgroundColor = styleRow(1:6, unlist(expected_colors[8:13, sample_headers[[1]]])),
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
        list(expected_data, expected_sample)
    )
    
    # also test without reference
    summary_stat_table <- summary_stat_table %>%
        select(-contains('reference'))
    sample_headers <- sample_headers[1]
    
    expected_data <- expected_tb[1:7,] %>%
        select(-3) %>%
        set_names(c('Data QC measures', sample_headers)) %>%
        datatable(
            options = list(
                dom = 't',
                pageLength = 7,
                rowCallback = JS(
                    "function(row, data, displayNum, displayIndex, dataIndex) {",
                    "let help_text = ", vector_to_js(summary_row_help[1:7]), ";",
                    "$('td', row).attr('title', help_text[data[0]]);",
                    "}"
                )
            ),
            rownames = FALSE
        ) %>% 
        # Color coding of values
        formatStyle(
            2, 
            backgroundColor = styleRow(1:7, unlist(expected_colors[1:7, sample_headers[[1]]])),
            textAlign = 'center'
        ) %>%
        formatStyle(
            3, 
            backgroundColor = styleRow(1:7, unlist(expected_colors[1:7, sample_headers[[1]]])),
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
        list(expected_data, expected_sample)
    )    
})

#SNV_table_output(SNV_table, roi_gr, snv_hotspot_tb, config, report_config, out_format = 'html', caption = NULL)
test_that("SNV_table_output", {
    
    test_config <- config
    test_config$settings[['SNV_analysis']] <- list(
        snv_hotspot_table = test_path('../data/minimal-snv-hotspots.tsv'),
        flag_GenCall_minimum = 0.2,
        variant_selection = list(
          Impact = list('HIGH', 'MODERATE'),
          Annotation_regex = NULL,
          include_all_ROI_overlaps = TRUE
        ),
        critical_SNV = list('ROI-overlap', 'hotspot-match'),
        reportable_SNV = list('hotspot-gene', 'protein-ablation', 'protein-changing'),
        protein_ablation_annotations = list(
            Impact = list('HIGH'),
            Annotation_regex = NULL
        ),
        protein_change_annotations = list(
            Impact = list(),
            Annotation_regex = 'missense_variant|inframe'
        )
    )
    report_config <- list()
    
    roi_gr <- tibble(
        list_name = 'ROI',
        hotspot = c('1:5000-5250', '1:8000-8250', 'DDX11L1'),
        mapping = c('position', 'position', 'gene_name'),
        call_type = 'any',
        check_score = 50,
        description = c('ROI_1: 1:5000-5250', 'ROI_2: 1:5000-5250', 'DDX11L1'),
        description_doi = NA_character_,
        description_htmllinks = c('ROI_1: 1:5000-5250', 'ROI_2: 1:5000-5250', 'DDX11L1'),
        order = 1:3,
    ) %>% 
        get_roi_gr(gr_genes, gr_info, target_chrom_style='UCSC') 
    
    snv_hotspot_tb <- load_hotspot_table(test_config, 'snv_hotspot')
    
    outcols <- c(
        'seqnames', 'start', 'REF', 'ALT', 'ID', 'FILTER', 'GT', 'GenTrain_Score', 'GenCall_Score',
        'Annotation', 'Impact', 'gene_name', 'Transcript_ID', 'Transcript_BioType',	'HGVS.c', 'HGVS.p', 'ROI_hits'
    )
    
    SNV_table <- sample_SNV_tb %>%
        dplyr::rename(Impact = Annotation_Impact) %>%
        select(all_of(outcols)) %>%
        filter(Impact %in% c('HIGH', 'MODERATE'))
    # Add some roi, Impact & annotation highlights
    SNV_table <- bind_rows(
        SNV_table %>% mutate(
            SNV_category = c('hotspot-gene', 'hotspot-gene', 'hotspot-match', 'hotspot-gene', 'hotspot-match'),
            SNV_label = c('reportable', 'reportable', 'critical', 'reportable', 'critical')
        ),
        SNV_table %>% mutate(
            SNV_label = c('unreliable impact', rep('reportable', 4)),
            SNV_category = c('ROI-overlap', 'protein-ablation', 'protein-changing', 'protein-ablation', 'protein-changing'),
            ROI_hits = c('DDX11L1', rep(NA, 4))
        )
    ) %>%        
        mutate(
            ref_GT = NA, # function will output ligcal NA vector
            ref_GenCall_Score = NA_real_,
            SNV_category = factor(SNV_category, levels = defined_labels$SNV_category_labels),
            SNV_label = factor(SNV_label, levels = defined_labels$SNV_labels),
        ) %>%
        arrange(SNV_label, SNV_category)
    
    # setup expected DT    
    column_help_text <- c(
        'Chromosome of the SNV/SNP',
        'Position of the SNV/SNP',
        'SNV/SNP in the format "Chr:Pos:REF>ALT"',
        'ID of the SNV/SNP from the illumina array.\\nNote: rsIDs may not always be reliable',
        paste0('Designation label for the SNV/SNP (', paste(defined_labels$SNV_labels, collapse = ', '), ')'),
        paste0('Evaluation category for the SNV/SNP (', paste(defined_labels$SNV_category_labels, collapse = ', '), ')'),
        'Reference allele of the SNV/SNP',
        'Alternative allele of the SNV/SNP',
        'Genotype of the SNV/SNP for 2 Allels: 0 stands for the reference allele, 1 for the alternative allele',
        'Reference sample genotype of the SNV/SNP for 2 Allels (0 - allele, 1 - alternative allele)',
        'Regions of interest overlapping with this SNV/SNP',
        'Gene name affected by the SNV/SNP',
        'Effect/Annotation of the SNV/SNP on the gene (from mehari)',
        'Impact annotation of the SNV/SNP on the gene (from mehari)',
        'Ensembl ID of the selected transcript for the affected gene',
        'Description of the selected transcript.\\nNote: ManeSelect designates the (medically) primary transcript of a gene.',
        'HGVS.c notation of the mutation/effect of the SNV/SNP on the transcript',
        'HGVS.p notation of the mutation/effect of the SNV/SNP on the protein',
        'GenTrain Score of the SNV/SNP',
        'GenCall Score of the SNV/SNP',
        'GenCall Score of the SNV/SNP in the reference sample'            
    )
    
    expected_dt <- SNV_table  %>%
        dplyr::rename(Chromosome = seqnames, Position = start) %>%
        mutate(
            desc = c('desc{1},{2}', 'desc', rep('SNV hotspot gene (see hotspot coverage)', 3), rep(NA, 5)),
            SNV = paste0(Chromosome, ':', format(Position, big.mark = '.', decimal.mark = ','), ':', REF, '>', ALT),
            ROI_hits = ifelse(1:dplyr::n() == 10, "<span class=\"badge badge-red\" title=\"ROI&#013;DDX11L1\">DDX11L1</span>", ROI_hits),
            HGVS.p = ifelse(
                1:dplyr::n() %in% 1:2, 
                paste0("<span class=\"badge badge-red\" title=\"test SNV hotspots&#013;", desc, '">', HGVS.p, "</span>"),
                HGVS.p
            ),
            gene_name = ifelse(
                1:dplyr::n() %in% 3:5, 
                paste0("<span class=\"badge badge-orange\" title=\"test SNV hotspots&#013;", desc, '">', gene_name, "</span>"),
                gene_name
            ),
            Impact = ifelse(
                1:dplyr::n() %in% 6:7, 
                paste0("<span class=\"badge badge-orange\">", Impact, "</span>"),
                Impact
            ), 
            Annotation = ifelse(
                1:dplyr::n() %in% 8:9, 
                paste0("<span class=\"badge badge-orange\">", Annotation, "</span>"),
                Annotation
            )
        ) %>%
        select(
            Chromosome, Position, SNV, ID, SNV_label, SNV_category, 
            REF, ALT, GT, ref_GT, ROI_hits, gene_name, Impact, Annotation,  
            Transcript_ID, Transcript_BioType, HGVS.c, HGVS.p,
            GenTrain_Score, GenCall_Score, ref_GenCall_Score
        ) %>%
        rename_with(format_column_names) %>%
        datatable(
            rownames = FALSE,
            escape = FALSE,
            extensions = c('Buttons', 'Scroller'),
            filter = 'top',
            caption = NULL,
            options = list(
                scrollY = 300, scrollCollapse = TRUE, scrollX =  TRUE, scroller = TRUE,
                dom = 'Bftilp',
                buttons = c('colvis', 'copy', 'csv', 'excel', 'print'),
                columnDefs = list(
                    list(targets = c(0:1,3,5:7,14:16,18,20), visible = FALSE)
                )
            ),
            callback = JS(
                "var info_text = ", vector_to_js(column_help_text), ";",
                "header = table.columns().header();",
                "for (var i = 0; i < info_text.length; i++) {",
                "  $(header[i]).attr('title', info_text[i]);",
                "};"
            )
        ) %>%
        formatRound(c('Position'), digits = 0, mark = '.')
    
    SNV_table_output(SNV_table, roi_gr, snv_hotspot_tb, test_config, report_config, defined_labels) %>%
        expect_equal(expected_dt)
    
    #TODO: add tests with empty roi_gr, empty snv_hotspot_tb & empty SNV_table
})