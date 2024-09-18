library(tidyverse)
library(plyranges)
library(DT)
library(knitr)
library(testthat)

source(test_path("../../stemcnv_check/scripts/R/helper_functions.R"))
source(test_path("../../stemcnv_check/scripts/R/R_table_functions.R"))

# Functions to test:
# - [ ] vector_to_js
# - [ ] format_column_names
# - [ ] simple_table_output
# - [ ] summary_table
# - [x] format_hotspots_to_badge
# - [ ] CNV_table_output
# - [ ] gene_table_output
# - [x] hotspot_table_output

# Test `format_hotspots_to_badge` function
# format_hotspots_to_badge <- function(hotspot_vec, CNVtype_vec, gene_details, listname = 'high_impact')
test_that("format_hotspots_to_badge", {
    testthat::local_edition(3)
    hotspot_vec <- c("", "12p13.3", "12p13.3", "TP53", "18q21|JAK2", "18q21|JAK2")
    CNVtype_vec <- c("gain", "gain", "LOH", "gain", "loss", "gain")
    # 1 - empty
    # 2 - gband hit, matching CNV (gain)
    # 3 - gband hit, not matching CNV (Note: this case should never occur; warning needed?)
    # 4 - gene hit
    # 5 - gene hit & gband hit matching CNV (loss)
    # 6 - gene hit & gband hit not matching CNV (Note: theoretically possible)

    gene_details <- read_tsv(
        test_path('../../stemcnv_check/supplemental-files/HighImpact-stemcell-hotspots.tsv'), 
        show_col_types = FALSE
    )
  
    expected <- c(
        '-', 
        '<span class="badge badge-HI" title="high_impact list name: StemCell-Hotspots&#013;Annotation source:&#013;StemCNV-check curation; ISCCR guidelines 2023&#013; doi.org/10.1038/ncomms5825 &#013; https://doi.org/10.1101/2021.05.22.445238">12p13.3</span>', 
        '12p13.3', 
        '<span class="badge badge-HI" title="high_impact list name: StemCell-Hotspots&#013;Annotation source:&#013;StemCNV-check curation; doi.org/10.1038/s41588-022-01147-3">TP53</span>', 
        '<span class="badge badge-HI" title="high_impact list name: StemCell-Hotspots&#013;Annotation source:&#013;StemCNV-check curation; ISCCR guidelines 2023">18q21</span><span class="badge badge-HI" title="high_impact list name: StemCell-Hotspots&#013;Annotation source:&#013;StemCNV-check curation; doi.org/10.1038/s41588-022-01147-3">JAK2</span>', 
        '18q21<span class="badge badge-HI" title="high_impact list name: StemCell-Hotspots&#013;Annotation source:&#013;StemCNV-check curation; doi.org/10.1038/s41588-022-01147-3">JAK2</span>'
    )
    expect_equal(
        format_hotspots_to_badge(hotspot_vec, CNVtype_vec, gene_details, 'high_impact'),
        expected
    )
    
    #test with include_hover = FALSE & listname = highlight
    expected <- c(
        '-', 
        '<span class="badge badge-HL">12p13.3</span>', 
        '12p13.3', 
        '<span class="badge badge-HL">TP53</span>', 
        '<span class="badge badge-HL">18q21</span><span class="badge badge-HL">JAK2</span>', 
        '18q21<span class="badge badge-HL">JAK2</span>'
    )
    expect_equal(
        format_hotspots_to_badge(hotspot_vec, CNVtype_vec, gene_details, 'highlight', FALSE),
        expected
    )
})

# hotspot_table_output(hotspots, plotsection, high_impact_tb, highlight_tb, report_config, out_format) %>%
test_that("hotspot_table_output", {
    hotspots <- c('TP53', '18q21')
    high_impact_tb <- read_tsv(
        test_path('../../stemcnv_check/supplemental-files/HighImpact-stemcell-hotspots.tsv'), 
        show_col_types = FALSE
    )
    highlight_tb <- tibble()
    # these aren't used so far
    plotsection <- 'test'
    report_config<- list()
    
    expected_tb <- high_impact_tb %>%
        filter(hotspot %in% hotspots) %>%
        select(hotspot, list_name, source, check_score, comment, mapping, call_type) %>%
        rename_with(format_column_names) 
    
    expected <- expected_tb %>%
        mutate(Hotspot = paste0('<span class="badge badge-HI">', Hotspot, '</span>')) %>%
        datatable(
            options = list(
                dom = 'Bt', 
                extensions = c('Buttons'),
                buttons = c('colvis', 'copy', 'print'),
                pageLength = 2,
                columnDefs = list(list(targets = 5:6, visible = FALSE))
            ),
            rownames = FALSE,
            escape = FALSE
        )
    
    hotspot_table_output(hotspots, plotsection, high_impact_tb, highlight_tb, report_config, 'html') %>%
        expect_equal(expected)
    
    # test non-html output
    expected <- expected_tb %>% 
        select(1:5) %>%
        kable()
    hotspot_table_output(hotspots, plotsection, high_impact_tb, highlight_tb, report_config, 'not-html') %>%
        expect_equal(expected)
    
    # test with highlight table
    highlight_tb <- high_impact_tb %>%
        filter(hotspot == '18q21') %>%
        mutate(list_name = 'highlight')
    high_impact_tb <- high_impact_tb %>%
        filter(hotspot != '18q21')
    expected <- expected_tb %>% 
        mutate(`List Name` = ifelse(Hotspot == '18q21', 'highlight', `List Name`)) %>%
        select(1:5) %>%
        kable()
    hotspot_table_output(hotspots, plotsection, high_impact_tb, highlight_tb, report_config, 'not-html') %>%
        expect_equal(expected)
})