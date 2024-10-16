library(tidyverse)
library(DT)
library(knitr)

vector_to_js <- function(v) {
    if (is.null(names(v))) {
        return(paste0("['", paste(v, collapse = "','"), "']"))
    } else {
        return(paste0('{', paste0(paste0("'", names(v), "': '", v, "'"), collapse = ', '), '}'))
    }
}

format_column_names <- function(n) {
    str_replace_all(n, '_', ' ') %>%
        str_to_title() %>%
        str_replace('Cnv', 'CNV') %>%
        str_replace('Loh', 'LOH') %>%
        str_replace('Snp', 'SNP') %>% 
        str_replace('Lrr', 'LRR') %>%
        str_replace('Roi', 'ROI') %>%
        str_replace('Id', 'ID')
}

simple_table_output <- function(tb, out_format='html', caption=NULL) {
    if (out_format == 'html') {
        return(
            datatable(
                tb,
                caption = caption,
                options = list(dom = 't', pageLength = nrow(tb)),
                rownames = FALSE
            )
        )
    } else {
        return(kable(tb, caption = caption))
    }
}

summary_table <- function(summary_stat_table, sample_headers, config) {
    
    Combined.metrics <- summary_stat_table %>%
        select(-ends_with('eval')) %>%
        filter(Description != 'sample_id') %>%
        mutate(Description = format_column_names(Description)) %>%
        set_names(c(' ', sample_headers))
    
    green_color <- ifelse(params$out_format == 'html', 'rgb(146,208,80)', 'green') #rgb(146,208,80) // lightgreen
    Combined.colors <- summary_stat_table %>%
        filter(Description != 'sample_id') %>%
        select(ends_with('eval')) %>%
        mutate(
            across(
                everything(),
                ~ case_when(
                    . == 'critical' ~ 'red',
                    . == 'warning' ~ 'orange',
                    . == 'unusual' ~ 'yellow',
                    . == 'OK' ~ green_color,
                    TRUE ~ 'white'
                )
            )
        ) %>%
        set_names(sample_headers)
        
    if (params$out_format == 'html') {
    
        call_count_excl_filters <- config$evaluation_settings$summary_stat_warning_levels$call_count_excl_filters
        ignored_calls <- ifelse(
            is.null(call_count_excl_filters) || length(call_count_excl_filters) == 0,
            '',
            paste0('\\nCalls with any of these Filters are not counted: ', call_count_excl_filters %>% paste(collapse = '|'))
        )

        summary_row_help <- c(
            "Call Rate" = 'The Illumina call rate corresponds to the percentage of probes for which a clear genotype could be determined.\\nLow values are strong indicator of sample quality issues that also impact make any further analysis including CNV calling.',
            "Computed Gender" = 'The Illumina computed gender (sex) based on X and Y chromosome probes.\\nMismatches with annotated sex can indicate annotation mistakes, sample swaps, or severe quality issues.',
            "SNPs Post Filter" = 'The percentage of SNP probes retained after the employed StemCNV-check filter strategy.',
            "SNP Distance To Reference" = 'The number of probes with a different genotype than the reference sample.\\nIncreased values can indicate a sample swap or a considerable number of genomic mutations between sample and references.',
            "Loss Gain Log2ratio" = paste0(
                'The log2 transformed ratio of loss and gain CNV calls.\\n',
                'Deviation from equal balance (0) can indicate potential quality issues or problems with CNV calling.',
                ignored_calls
            ),
            "Total Calls CNV" = paste0('The total number of CNV (gain/loss) calls.', ignored_calls),
            "Total Calls LOH" = paste0('The total number of LOH calls.', ignored_calls)
        )
        if ("Reportable Calls CNV" %in% Combined.metrics$` `) {
            summary_row_help <- c(summary_row_help,
                "Reportable Calls CNV" = 'The number of CNV calls designated as "reportable".',
                "Reportable Calls LOH" = 'The number of LOH calls designated as "reportable".'
            )
        }
        if ("Critical Calls CNV" %in% Combined.metrics$` `) {
            summary_row_help <- c(summary_row_help,
                "Critical Calls CNV" = 'The number of CNV calls designated as "critical".',
                "Critical Calls LOH" = 'The number of LOH calls designated as "critical".'
            )
        }
        datatable(Combined.metrics,
            options = list(
                dom = 't',
                pageLength = nrow(Combined.metrics),
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
            backgroundColor = styleRow(1:nrow(Combined.metrics), unlist(Combined.colors[, sample_headers[[1]]])),
            textAlign = 'center'
        ) %>%
        # DT can take specification of non-existing columns - just need a workaround for the sample_header call
        formatStyle(
            3,
            backgroundColor = styleRow(1:nrow(Combined.metrics), unlist(Combined.colors[, sample_headers[[length(sample_headers)]]])), 
            textAlign = 'center'
        ) %>%
        # Make all critical (= potentially red) rows have bold text
        formatStyle(
            1, 
            fontWeight = styleEqual(
                unlist(config$evaluation_settings$summary_stat_warning_levels$last_level_critical) %>% format_column_names(), 
                'bold'
            )
        )
        
    } else {
    
        tbout <- kable(Combined.metrics, align = c('l', rep('c', length(sample_headers))), format = 'latex') %>%
            column_spec(2, background = unlist(Combined.colors[, sample_headers[[1]]]))
        if (!is.na(ref_id)) {
           tbout <- column_spec(tbout, 3, background = unlist(Combined.colors[, sample_headers[[2]]]))
        }
        tbout
    }


}

format_hotspots_to_badge <- function(
    hotspot_vec, CNVtype_vec, hotspot_table, listname = 'high_impact', include_hover = TRUE
) {
    if (listname == "high_impact") {
        shorthand <- 'HI'
    } else if (listname == "highlight") {
        shorthand <- 'HL'
    } else {
        stop(str_glue("Unsupported list type '{listname}', only 'high_impact' and 'highlight' are defined"))
    }

    hotspot_table.any <- hotspot_table %>%
        mutate(call_type = str_replace(call_type, 'any', 'any;gain;loss;LOH')) %>%
        separate_rows(call_type, sep = ';')
    
    tb <- tibble(
        hotspot = str_split(hotspot_vec, '\\|'),
        call_type = CNVtype_vec
    ) %>%
        mutate(id = 1:dplyr::n()) %>%
        unnest(hotspot) %>%
        # simple joining with base hotspot table by call_type does not work due to 'any' not matching
        left_join(hotspot_table.any, by = c('hotspot', 'call_type')) %>%
        rowwise() %>%
        mutate(
            include_hover = include_hover,
            description = str_replace_all(description, '\\\\n', '&#013;') %>%
                str_replace_all('\\n', '&#013;'),
            out_str = ifelse(
                hotspot != "" & !is.na(list_name),
                paste0(
                    str_glue('<span class="badge badge-{shorthand}"'),
                    ifelse(
                        include_hover,
                        paste0(
                            str_glue(' title="{list_name}&#013;'),
                            ifelse(!is.na(check_score), str_glue('Check_Score contribution: {check_score}&#013;'), ''),
                            description,
                            '"'
                        ),
                        ''
                    ),                    
                    str_glue('>{hotspot}</span>')
                ),
                hotspot
            )
        ) %>%
        group_by(id) %>%
        summarise(
            sep = ifelse(all(str_detect(out_str, '^<')), "", ', '),
            hotspot = base::paste(out_str, collapse='')
        ) %>%
        mutate(hotspot = ifelse(hotspot %in% c('', 'NA'), '-', hotspot))

    if (nrow(tb) == 0 | all(is.na(hotspot_vec) | hotspot_vec == 'NA')) {
        out_vec <- ifelse(is.na(hotspot_vec), '-', hotspot_vec)
        return(out_vec)
    } else {
        return(tb$hotspot)
    }
}

CNV_table_output <- function(tb, plotsection, high_impact_tb, highlight_tb, gr_info, report_config, caption = NULL) {
    always_include <- report_config$call.data.and.plots[[plotsection]]$always_include
    # Reorder & subset columns
    tb <- tb %>%
    # The links for internal plots do switch to the correct plot, but switching the active marked tab requires more JS, sadly the Rmd tabs do not have IDs making this very tricky
        mutate(
            Plot = ifelse(
                i <= report_config$call.data.and.plots[[plotsection]]$min_number_plots | Call_label %in% always_include,
                paste0(
                    '<a data-toggle="tab" href="#',
                    pmap_chr(., CNV_ID_str),
                    '">Nr. ', i, '</a>'
                ),
                paste0(
                    '<a href="" onclick="window.open(\'./',
                    image_folder, "/CNV_call.", plotsection, ".nr", i, '.plot-1.png\'); return false;">',
                    'Nr. ', i,' (ext. plot)</a>'
                )
            ), 
            Precision_Estimate = ifelse(is.na(Precision_Estimate), '-', as.character(Precision_Estimate)),
            HighImpact = map2_chr(HighImpact, CNV_type, \(hi,c) format_hotspots_to_badge(hi, c, high_impact_tb, 'high_impact')),
            Highlight = map2_chr(Highlight, CNV_type, \(hi,c) format_hotspots_to_badge(hi, c, highlight_tb, 'highlight')),
            genome_bands = pmap_chr(., \(chrom, start, end, ...) {
                gr_info %>% 
                    filter_by_overlaps(GRanges(seqnames = chrom, strand = '*', ranges = IRanges(start, end))) %>%
                    as_tibble() %>%
                    summarise(
                        gbands = ifelse(
                            dplyr::n() > 1,
                            paste(section_name[start == min(start)], '-', section_name[start == max(start)]),
                            section_name
                        )
                    ) %>%
                    pull(gbands)
            })
        ) %>%
        dplyr::rename(copynumber = CN) %>%  
        select(
            sample_id, ID, i, #invis 0-2
            Plot, Call_label, Check_Score,
            CNV_type, chrom, Size, genome_bands,
            start, end, #invis 10-11
            CNV_caller, HighImpact, Highlight, ROI_hits,
            Precision_Estimate, probe_coverage_gap, high_probe_density,
            # invis: 19++
            copynumber, LRR, n_probes, n_uniq_probes, #n_premerged_calls, caller_confidence,
            caller_merging_coverage, Gap_percent
        ) 

    if (params$out_format == 'html') {
      
        #TODO adjust plot/link hovertext to match settings  
        column_help_text <- c(
            'Sample_ID from the input "sample_table"',
            'Internal ID for the CNV call',
            'Number of the CNV call, sorted by descending Check_Score',
            'Link to the plot of the CNV call\\nNote: For the Top20/critical CNVs clicking on the link will switch the active plot below. For other CNVs it will open the plot in a new browser tab.',
            'Designation label for the CNV call, (Critical, Reportable, Reference genotype, ROI)',
            'Check_Score of the CNV call, calculated based on size, overlap high impact or highlight list, or other genes',
            'Type of CNV call (gain, loss, LOH)',
            'Chromosome of the CNV call',
            'Size of the CNV call (in base pairs)',
            'Genome bands overlapping the CNV call',
            'Start position of the CNV call',
            'End position of the CNV call',
            'Caller tools detecting this CNV call',
            'High impact hotspots overlapping with this CNV call',
            'Highlight hotspots overlapping with this CNV call',
            'Regions of interest overlapping with this CNV call',
            #FIXME (future): add a doi for precision benchmark once available
            'Precision estimate of the CNV call, based on internal benchmarking',
            'Call has a gap in probe coverage\\nBased on percentage of call without probes and size of the call. See config min.perc.gap_area and gap_area.uniq_probes.rel for details',
            'Call has higher probe density than {density.quantile.cutoff} percent of the the array',
            '(Estimated) copy number of the CNV call',
            'Median Log R Ratio of the CNV call',
            'Number of SNP probes (post filtering) in the CNV call area',
            'Number of unique probe positions (post filtering) in the CNV call area',
            # 'Number of calls from the caller before merging, comma separated for mutliple callers',
            # 'Call confidence of the caller where available, comma separated for multiple callers',
            'For calls with multiple CNV callers: percentage overlap of each caller with the combined call area',
            'Percentage of the CNV call area with a gap probe coverage'
        )

        dt <- datatable(
            tb %>% rename_with(format_column_names),
            rownames = FALSE,
            escape = FALSE,
            extensions = c('Buttons', 'Scroller'),
            filter = 'top',
            caption = caption,
            options = list(
                scrollY = 300, scrollCollapse = TRUE, scrollX =  TRUE, scroller = TRUE,
                dom = 'Bftilp',
                buttons = c('colvis', 'copy', 'csv', 'excel', 'print'),
                columnDefs = list(
                    #This uses 0-indexing vs the usual R 1-indexing
                    list(targets = c(0:2,10:11,19:(ncol(tb)-1)), visible = FALSE)
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
        formatRound(c('Start', 'End', 'Size'), digits = 0, mark = '.') %>%
        formatRound(c('Check Score', 'Gap Percent'), digits = 2)
    return(dt)
    } else {
        tb <- tb %>% 
            select(
                CNV_type, Check_Score,
                chrom, start, end, Size, genome_bands, CNV_caller,
                HighImpact, Highlight,
                Precision_Estimate, probe_coverage_gap, high_probe_density
            ) %>% 
            rename_with(format_column_names)
        return(kable(tb, caption = caption))
    }
}

gene_table_output <- function(
    tb, plotsection, high_impact_tb, highlight_tb, report_config, caption = NULL, extra_cols = c()
) {

    if (report_config$call.data.and.plots[[plotsection]]$include.gene.table.details == 'Call') {
        tb <- filter(tb, direct_hit)
        gene_area <- 'call'
    } else {
        extra_cols <- c(extra_cols, 'direct_hit')
        gene_area <- 'plot'
    }
    if (nrow(tb) == 0) {
        return(str_glue('No genes in the {gene_area} area.'))
    }
    tb <- tb %>%
        mutate(
            across(any_of(c('direct_hit', 'gene_type', 'strand')), ~ factor(.)),
            high_impact = factor(ifelse(high_impact, 'hit', '-'), levels = c('hit', '-')),
            highlight = factor(ifelse(highlight, 'hit', '-'), levels = c('hit', '-')),
            name_is_geneid = str_detect(gene_name, 'ENSG[0-9]{11}'),
            # REEV: gene_id *should* work, but won't if they are deprectated/not in Annonars
            REEV = str_glue("<a href='https://reev.bihealth.org/gene/{gene_name}' target='_blank' rel='noopener noreferrer'>{gene_name}</a>"),
            GTex = str_glue("<a href='https://gtexportal.org/home/gene/{gene_name}' target='_blank' rel='noopener noreferrer'>{gene_name}</a>"),
            NCBI = ifelse(name_is_geneid, '-', str_glue("<a href='https://pubmed.ncbi.nlm.nih.gov/?term={gene_name}' target='_blank' rel='noopener noreferrer'>{gene_name}</a>")),
            Ensembl = str_glue("<a href=' https://www.ensembl.org/Homo_sapiens/Gene/Summary?g={gene_id}' target='_blank' rel='noopener noreferrer'>{gene_id}</a>"),
            #Reformat gene name
            gene_name = ifelse(
                high_impact == 'hit',
                map2_chr(gene_name, CNVtype, \(g, c) format_hotspots_to_badge(g,c, high_impact_tb,'high_impact', FALSE)),
                map2_chr(gene_name, CNVtype, \(g, c) format_hotspots_to_badge(g,c, highlight_tb, 'highlight', FALSE))
            ),
        ) %>%
        arrange(high_impact, highlight, desc(direct_hit), start) %>%
        select(
            gene_name, gene_id, seqnames, start, end, strand, high_impact, highlight,
            any_of(extra_cols), REEV, GTex, NCBI, Ensembl
        )
    if (params$out_format == 'html') {
        colors1 <- ifelse(tb$high_impact == 'hit', 'red' , 'white')
        colors2 <- ifelse(tb$highlight == 'hit', 'orange', 'white')
        dt <- datatable(
            tb, rownames = FALSE, escape = FALSE,
            options = list(
                dom = 'Bftilp', pageLength = 10,
                extensions = c('Buttons'),
                buttons = c('colvis', 'copy', 'csv', 'excel', 'print'),
                columnDefs = list(list(targets = 1:7, visible = FALSE))
                )
            ) %>%
            formatStyle('high_impact', backgroundColor = styleRow(1:nrow(tb), colors1)) %>% # textAlign = 'center'
            formatStyle('highlight', backgroundColor = styleRow(1:nrow(tb), colors2))
        return(dt)
    } else {
        tb <- tb %>% select(gene_name, gene_id, high_impact, highlight, any_of(extra_cols))
        return(kable(tb, caption = caption))
    }
}

hotspot_table_output <- function(
    hotspots, cnv_type, plotsection, high_impact_tb, highlight_tb, report_config, out_format, caption = NULL
){
    tb <- bind_rows(
        high_impact_tb,
        highlight_tb
    ) %>%
        filter(hotspot %in% hotspots & call_type %in% c('any', cnv_type))
    
    if (out_format == 'html') {
        dt <- datatable(
            tb %>%
                select(-description) %>%
                dplyr::rename(description = description_htmllinks) %>%
                select(hotspot, call_type, list_name, description, check_score, any_of(colnames(tb))) %>%
                mutate(
                    hotspot = ifelse(
                        list_name %in% unique(high_impact_tb$list_name),
                        map2_chr(hotspot, call_type, \(g, c) format_hotspots_to_badge(g,c, high_impact_tb,'high_impact', FALSE)),
                        map2_chr(hotspot, call_type, \(g, c) format_hotspots_to_badge(g,c, highlight_tb, 'highlight', FALSE))
                    ),
                    description = str_replace_all(description, '&#013;', '<br/>')
                ) %>%
                rename_with(format_column_names),
            caption = caption,
            options = list(
                dom = 'Bt', 
                extensions = c('Buttons'),
                buttons = c('colvis', 'copy', 'print'),
                pageLength = nrow(tb),
                # DT cols are 0-indexed, 1 col removed
                columnDefs = list(list(targets = 5:(ncol(tb)-2), visible = FALSE))
            ),
            rownames = FALSE,
            escape = FALSE
        )
        return(dt)
    } else {
        return(
            kable(
                tb %>%
                    dplyr::rename(dois = description_doi) %>%
                    select(hotspot, call_type, list_name, description, check_score, dois) %>%
                    rename_with(format_column_names),
                caption = caption
            )
        )
    }
}