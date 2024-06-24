
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
    str_replace('Roi', 'ROI') %>%
    str_replace('Id', 'ID')
}

simple_table_output <- function(tb, caption=NULL) {
  if (params$out_format == 'html') {
    return(datatable(tb, caption = caption,
                     options = list(dom = 't', pageLength = nrow(tb)), rownames = FALSE))
  } else {
    return(kable(tb, caption = caption))
  }
}

summary_table <- function(Combined.metrics, Combined.colors, sample_headers) {
  if (params$out_format == 'html') {

    summary_row_help <- c(
            "Call Rate" = 'The Illumina call rate corresponds to the percentage of probes for which a clear genotype could be determined.\\nLow values are strong indicator of sample quality issues that also impact make any further analysis including CNV calling.',
            "Computed Gender" = 'The Illumina computed gender (sex) based on X and Y chromosome probes.\\nMismatches with annotated sex can indicate annotation mistakes, sample swaps, or severe quality issues.',
            "SNPs Post Filter" = 'The percentage of SNP probes retained after the employed StemCNV-check filter strategy.',
            "SNP Distance To Reference" = 'The number of probes with a different genotype than the reference sample.\\nIncreased values can indicate a sample swap or a considerable number of genomic mutations between sample and references.',
            "Loss Gain Log2ratio" = 'The log2 transformed ratio of loss and gain CNV calls.\\nDeviation from equal balance (0) can indicate potential quality issues or problems with CNV calling.',
            "Total Calls CNV" = 'The total number of CNV (gain/loss) calls.',
            "Total Calls LOH" = 'The total number of LOH calls.'
    )
    #TODO: add actual tresholds into the help text
    if (!str_detect(not_enabled_label, 'reportable')) {
      summary_row_help <- c(summary_row_help,
             "Reportable Calls CNV" = 'The number of CNV calls with an Impact Score above the reportable threshold.',
             "Reportable Calls LOH" = 'The number of LOH calls with an Impact Score above the reportable threshold.'
      )
    }
    if (!str_detect(not_enabled_label, 'critical')) {
      summary_row_help <- c(summary_row_help,
             "Critical Calls CNV" = 'The number of CNV calls with an Impact Score above the critical threshold.',
             "Critical Calls LOH" = 'The number of LOH calls with an Impact Score above the critical threshold.'
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
                  "}")
              ),
              rownames = FALSE) %>%
        formatStyle(2, backgroundColor = styleRow(1:nrow(Combined.metrics), unlist(Combined.colors[, sample_headers[[1]]])), textAlign = 'center') %>%
        # DT can take specification of non-existing columns - just need a workaround for the sample_header call
        formatStyle(3, backgroundColor = styleRow(1:nrow(Combined.metrics), unlist(Combined.colors[, sample_headers[[length(sample_headers)]]])), textAlign = 'center')
  } else {

    tbout <- kable(Combined.metrics, align = c('l', rep('c', length(sample_headers))), format = 'latex') %>%
      column_spec(2, background = unlist(Combined.colors[, sample_headers[[1]]]))
    if (!is.na(ref_id)) {
       tbout <- column_spec(tbout, 3, background = unlist(Combined.colors[, sample_headers[[2]]]))
    }
    tbout
  }


}

format_hotspots_to_badge <- function(hotspot_vec, CNVtype_vec, gene_details, listname = 'high_impact') {
  if (listname == "high_impact") {
    shorthand <- 'HI'
  } else if (listname == "highlight") {
    shorthand <- 'HL'
  } else {
    stop(str_glue("Unsupported list type '{listname}', only 'high_impact' and 'highlight' are defined"))
  }

  tb <- tibble(hotspot = str_split(hotspot_vec, ','),
               CNV_type = CNVtype_vec) %>%
          mutate(id = 1:dplyr::n()) %>%
          unnest(hotspot) %>%
          left_join(gene_details) %>%
          mutate(
            do_format = case_when(
              is.na(call_type) ~ FALSE,
              call_type == CNV_type ~ TRUE,
              call_type == 'any' ~ TRUE,
              TRUE ~ FALSE
            ),
            source = str_replace_all(source, '\\\\n', '&#013;') %>%
              str_replace_all('\\n', '&#013;'),
            out_str = ifelse(
              hotspot != "" & !is.na(list_name) & do_format,
              paste0(str_glue('<span class="badge badge-{shorthand}" title="'),
                     str_glue('{listname} list name: {list_name}&#013;'), #\\n
                     ifelse(!is.na(impact_score), str_glue('custom impact_score: {impact_score}&#013;'), ''),
                     str_glue('Annotation source:&#013;{source}'),
                     str_glue('">{hotspot}</span>')),
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

CNV_table_output <- function(tb, plotsection, high_impact_tb, highlight_tb, caption = NULL) {
  always_include <- report.setting('call.data.and.plots', plotsection, 'always_include')
  # Reorder & subset columns
  tb <- tb %>%
    # The links for internal plots do switch to the correct plot, but switching the active marked tab requires more JS, sadly the Rmd tabs do not have IDs making this very tricky
    mutate(Plot = ifelse(
             i <= report.setting('call.data.and.plots', plotsection, 'min_number_plots') | Call_Label %in% always_include,
             paste0('<a data-toggle="tab" href="#', pmap_chr(., CNV_ID_str, formatting = 'link'),
             #       ' onclick="$(', "'#", pmap_chr(., CNV_ID_str, formatting = 'link'), "').trigger('click')",
                    '">Nr. ', i, '</a>'),
             #       '">', pmap_chr(., CNV_ID_str), '</a>'),
             paste0('<a href="" onclick="window.open(\'./', image_folder, "/CNV_call.", plotsection, ".nr", i, '.plot-1.png\'); return false;">',
                    'Nr. ', i,' (ext. plot)</a>')
              # pmap_chr(., CNV_ID_str),' (ext. plot)</a>')
             ),
           Precision_Estimate = ifelse(is.na(Precision_Estimate), '-', as.character(Precision_Estimate)),
           high_impact_hits = map2_chr(high_impact_hits, CNV_type, \(hi,c) format_hotspots_to_badge(hi, c, high_impact_tb, 'high_impact')),
           highlight_hits = map2_chr(highlight_hits, CNV_type, \(hi,c) format_hotspots_to_badge(hi, c, highlight_tb, 'highlight')),
    ) %>%
    select(sample_id, ID, i, #invis 0-2
           Plot, Call_Label, Impact_Score,
           CNV_type, Chr, Size,
           Start, End, #invis 9-10
           CNV_caller, high_impact_hits, highlight_hits, ROI_hits,
           Precision_Estimate, probe_coverage_gap, high_probe_density,
           # invis: 18++
           copynumber, n_snp_probes, n_uniq_probe_positions, n_premerged_calls, caller_confidence,
           caller_merging_coverage,
           percent_gap_coverage )

  if (params$out_format == 'html') {

    column_help_text <- c(
            'Sample_ID from the input "sample_table"',
            'Internal ID for the CNV call',
            'Number of the CNV call, sorted by descending Impact Score',
            'Link to the plot of the CNV call\\nNote: For the Top20/critical CNVs clicking on the link will switch to the active plot below, without hightlighting the correct tab. For other CNVs it will open the plot in a new browser tab.',
            'Designation label for the CNV call, (Critical, Reportable, Reference genotype, ROI)',
            'Impact score of the CNV call, calculated based on size, overlap high impact or highlight list, or other genes',
            'Type of CNV call (gain, loss, LOH)',
            'Chromosome of the CNV call',
            'Size of the CNV call (in base pairs)',
            'Start position of the CNV call',
            'End position of the CNV call',
            'Caller tools detecting this CNV call',
            'High impact hotspots overlapping with this CNV call',
            'Highlight hotspots overlapping with this CNV call',
            'Regions of interest overlapping with this CNV call',
            #TODO: add a doi for precision benchmark once available
            'Precision estimate of the CNV call, based on internal benchmarking',
            'Call has a gap in probe coverage\\nBased on percentage of call without probes and size of the call. See config min.perc.gap_area and gap_area.uniq_probes.rel for details',
            'Call has higher probe density than {density.quantile.cutoff} percent of the the array',
            '(Estimated) copy number of the CNV call',
            'Number of SNP probes (post filtering) in the CNV call area',
            'Number of unique probe positions (post filtering) in the CNV call area',
            'Number of calls from the caller before merging, comma separated for mutliple callers',
            'Call confidence of the caller where available, comma separated for multiple callers',
            'For calls with multiple CNV callers: percentage overlap of each caller with the combined call area',
            'Percentage of the CNV call area with a gap probe coverage'
    )


    dt <- datatable(tb %>% rename_with(format_column_names),
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
                list(targets = c(0:2,9:10,18:(ncol(tb)-1)), visible = FALSE)
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
      formatRound(c('Impact Score'), digits = 2, mark='.')
    return(dt)
  } else {
    tb <- tb %>% select(CNV_type, Impact_Score,
                        Chr, Start, End, Size, CNV_caller,
                        high_impact_hits, highlight_hits,
                        Precision_Estimate, probe_coverage_gap, high_probe_density
    ) %>% rename_with(format_column_names)
    return(kable(tb, caption = caption))
  }
}

gene_table_output <- function(tb, plotsection, high_impact_tb, highlight_tb, caption = NULL, extra_cols = c()) {

  if (report.setting('call.data.and.plots', plotsection, 'include.gene.table.details') == 'Call') {
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
              map2_chr(gene_name, CNVtype, \(g, c) format_hotspots_to_badge(g,c, high_impact_tb,'high_impact')),
              map2_chr(gene_name, CNVtype, \(g, c) format_hotspots_to_badge(g,c, highlight_tb, 'highlight'))
              ),
    ) %>%
    arrange(high_impact, highlight, desc(direct_hit), start) %>%
    select(gene_name, gene_id, seqnames, start, end, high_impact, highlight, any_of(extra_cols), REEV, GTex, NCBI, Ensembl)
  if (params$out_format == 'html') {
    colors1 <- ifelse(tb$high_impact == 'hit', 'red' , 'white')
    colors2 <- ifelse(tb$highlight == 'hit', 'orange', 'white')
    dt <- datatable(tb, rownames = FALSE, escape = FALSE,
                     options = list(dom = 'Bftilp', pageLength = 10,
                                    extensions = c('Buttons'),
                                    buttons = c('colvis', 'copy', 'csv', 'excel', 'print'),
                                    columnDefs = list(list(
                                      targets = 1:6, visible = FALSE
                                    ))
                     )) %>%
        formatStyle('high_impact', backgroundColor = styleRow(1:nrow(tb), colors1)) %>% # textAlign = 'center'
        formatStyle('highlight', backgroundColor = styleRow(1:nrow(tb), colors2))
    return(dt)
  } else {
    tb <- tb %>% select(gene_name, gene_id, high_impact, highlight, any_of(extra_cols))
    return(kable(tb, caption = caption))
  }
}
