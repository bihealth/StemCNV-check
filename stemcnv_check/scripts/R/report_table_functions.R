library(tidyverse)
library(DT)
library(knitr)

# Helper function, could be defined elsewhere?
get_defined_labels <- function(config) {
    label_def_file <- file.path(config$snakedir, 'control_files', 'label_name_definitions.yaml')
    return(yaml.load_file(label_def_file))
}

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
        str_replace('Snv', 'SNV') %>%
        str_replace('Gt', 'GT') %>% 
        str_replace('Cn', 'CN') %>% 
        str_replace('Cdna', 'cDNA') %>% 
        str_replace('Cds', 'CDS') %>% 
        str_replace('Aa', 'AA') %>% 
        str_replace('Hgvs', 'HGVS') %>%
        str_replace('Lrr', 'LRR') %>%
        str_replace('Roi', 'ROI') %>%
        str_replace('Id', 'ID')
}

simple_table_output <- function(tb, out_format='html', caption=NULL, escape=TRUE) {
    if (out_format == 'html') {
        return(
            datatable(
                tb,
                caption = caption,
                escape = escape,
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
    
    sample_levels <- get_defined_labels(config)$sample_levels
    
    # R/html default green is very neon, use a nicer shade
    green_replace <- ifelse(params$out_format == 'html', 'rgb(146,208,80)', 'green') #rgb(146,208,80) // lightgreen
    Combined.colors <- summary_stat_table %>%
        filter(Description != 'sample_id') %>%
        select(ends_with('eval')) %>%
        mutate(
            across(
                everything(),
                ~ ifelse(. %in% names(sample_levels), sample_levels[[.]], 'white') %>%
                    str_replace('green', green_replace)
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
        summary_row_help <- c(summary_row_help,
            "Critical SNVs" = 'The number of detected SNVs designated as "critical".'
        )
        
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
        # Make all last level (= potentially red) rows have bold text
        formatStyle(
            1, 
            fontWeight = styleEqual(
                unlist(config$evaluation_settings$summary_stat_warning_levels$use_last_level) %>% format_column_names(), 
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
    hotspot_vec, CNVtype_vec, hotspot_table, color = 'orange', include_hover = TRUE
) {
    
    if (color %!in% c('orange', 'red')) {
        stop(str_glue("Unsupported color '{color}', only 'red' and 'orange' are defined"))
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
                    str_glue('<span class="badge badge-{color}"'),
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

CNV_table_output <- function(
    tb, plotsection, stemcell_hotspot_tb, dosage_sensitive_gene_tb, cancer_gene_tb, gr_info, report_config, caption = NULL
) {
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
            stemcell_hotspot = map2_chr(stemcell_hotspot, CNV_type, \(hi,c) format_hotspots_to_badge(hi, c, stemcell_hotspot_tb, 'red')),
            dosage_sensitive_gene = map2_chr(dosage_sensitive_gene, CNV_type, \(hi,c) format_hotspots_to_badge(hi, c, dosage_sensitive_gene_tb, 'orange')),
            cancer_gene = map2_chr(cancer_gene, CNV_type, \(hi,c) format_hotspots_to_badge(hi, c, cancer_gene_tb, 'orange')),
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
            CNV_caller, stemcell_hotspot, dosage_sensitive_gene, cancer_gene, ROI_hits,
            Precision_Estimate, probe_coverage_gap, high_probe_density,
            # invis: 20++
            copynumber, LRR, n_probes, n_uniq_probes, #n_premerged_calls, caller_confidence,
            caller_merging_coverage, Gap_percent
        ) 

    if (params$out_format == 'html') {
        
        column_help_text <- c(
            'Sample_ID from the input "sample_table"',
            'Internal ID for the CNV call',
            'Number of the CNV call, sorted by descending Check_Score',
            paste0(
                'Link to the plot of the CNV call\\nNote: For the Top',
                report_config$call.data.and.plots[[plotsection]]$min_number_plots, 
                ' or  all ',
                paste( 
                    always_include,
                    collapse = '/'
                ),
                ' CNVs clicking on the link will switch the active plot below. ',
                'For other CNVs it will open an externally saved plot in a new browser tab.'
            ),
            paste0(
                'Designation label for the CNV call, (',
                paste(names(config$evaluation_settings$CNV_call_labels), collapse = ', '),
                ')'
            ),
            'Check_Score of the CNV call, calculated based on size, overlap with stemcell_hotspot, dosage_sensivity or cancer_gene list, or other genes',
            'Type of CNV call (gain, loss, LOH)',
            'Chromosome of the CNV call',
            'Size of the CNV call (in base pairs)',
            'Genome bands overlapping the CNV call',
            'Start position of the CNV call',
            'End position of the CNV call',
            'CNV caller tools detecting this CNV call',
            'Stemcell hotspots overlapping with this CNV call',
            'Dosage sensitive genes overlapping with this CNV call',
            'Cancer genes overlapping with this CNV call',
            'Regions of interest overlapping with this CNV call',
            #FIXME (future): add a doi for precision benchmark once available
            'Precision estimate of the CNV call, based on internal benchmarking',
            paste0(
                'Call has a gap in probe coverage.\\nBased on percentage of call without probes and size of the call. ',
                'Based on `min.perc.gap_area` and `gap_area.uniq_probes.rel` from config settings:CNV_processing:call_processing'
            ),
            str_glue(
                'Call has higher probe density than {100*config$settings$CNV_processing$call_processing$density.quantile.cutoff}',
                ' percent of the the array (from config settings:CNV_processing:call_processing)'
            ),
            '(Estimated) copy number of the CNV call',
            'Median Log R Ratio of the CNV call',
            'Number of SNP probes (post filtering) in the CNV call area',
            'Number of unique probe positions (post filtering) in the CNV call area',
            # 'Number of calls from the caller before merging, comma separated for mutliple callers',
            # 'Call confidence of the caller where available, comma separated for multiple callers',
            'For calls with multiple CNV callers: percentage overlap of each caller with the combined call area',
            'Percentage of the CNV call area with a gap probe coverage'
        )
        # 
        # paste0(
        #     '<script type="text/javascript">',
        #     
        #     
        #     sep = '\n\n'
        #     
        # )

        dt <- datatable(
            tb %>% rename_with(format_column_names),
            rownames = FALSE,
            escape = FALSE,
            extensions = c('Buttons', 'Scroller', 'FixedColumns', 'Select'),
            filter = 'top',
            caption = caption,
            options = list(
                scrollY = 300, scrollCollapse = TRUE, scrollX =  TRUE, scroller = TRUE,
                dom = 'Bftilp',
                buttons = c('colvis', 'copy', 'csv', 'excel', 'print'),
                fixedColumns = list(leftColumns = 2),
                select = list(style = 'single', items = 'row', selector = 'td:first-child', togglebale = FALSE, rows = ''),
                columnDefs = list(
                    #This uses 0-indexing vs the usual R 1-indexing
                    list(targets = c(0:2,10:11,20:(ncol(tb)-1)), visible = FALSE)
                )
            ),
            callback = JS(
                "var info_text = ", vector_to_js(column_help_text), ";",
                "header = table.columns().header();",
                "for (var i = 0; i < info_text.length; i++) {",
                "  $(header[i]).attr('title', info_text[i]);",
                "};"
            ),
            selection = 'none'
        ) %>%
        formatRound(c('Start', 'End', 'Size'), digits = 0, mark = '.') %>%
        formatRound(c('Check Score', 'Gap Percent'), digits = 2)
    return(dt)
    } else {
        tb <- tb %>% 
            select(
                CNV_type, Check_Score,
                chrom, start, end, Size, genome_bands, CNV_caller,
                stemcell_hotspot, dosage_sensitive_gene, cancer_gene,
                Precision_Estimate, probe_coverage_gap, high_probe_density
            ) %>% 
            rename_with(format_column_names)
        return(kable(tb, caption = caption))
    }
}

gene_table_output <- function(
    tb, plotsection, stemcell_hotspot_tb, dosage_sensitive_gene_tb, cancer_gene_tb, report_config, caption = NULL, extra_cols = c()
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
            across(
                any_of(c('stemcell_hotspot', 'dosage_sensitive_gene', 'cancer_gene')),
                ~ factor(ifelse(., 'yes', '-'), levels = c('yes', '-'))
            ),
            # stemcell_hotspot = factor(ifelse(stemcell_hotspot, 'hit', '-'), levels = c('hit', '-')),
            # cancer_gene = factor(ifelse(cancer_gene, 'hit', '-'), levels = c('hit', '-')),
            name_is_geneid = str_detect(gene_name, 'ENSG[0-9]{11}'),
            # REEV: gene_id *should* work, but won't if they are deprectated/not in Annonars
            REEV = str_glue("<a href='https://reev.bihealth.org/gene/{gene_name}' target='_blank' rel='noopener noreferrer'>{gene_name}</a>"),
            GTex = str_glue("<a href='https://gtexportal.org/home/gene/{gene_name}' target='_blank' rel='noopener noreferrer'>{gene_name}</a>"),
            NCBI = ifelse(name_is_geneid, '-', str_glue("<a href='https://pubmed.ncbi.nlm.nih.gov/?term={gene_name}' target='_blank' rel='noopener noreferrer'>{gene_name}</a>")),
            Ensembl = str_glue("<a href=' https://www.ensembl.org/Homo_sapiens/Gene/Summary?g={gene_id}' target='_blank' rel='noopener noreferrer'>{gene_id}</a>"),
            #Reformat gene name
            gene_name = case_when(
                stemcell_hotspot == 'yes' ~ map2_chr(gene_name, CNVtype, \(g, c) format_hotspots_to_badge(g,c, stemcell_hotspot_tb,'red', FALSE)),
                dosage_sensitive_gene == 'yes' ~ map2_chr(gene_name, CNVtype, \(g, c) format_hotspots_to_badge(g,c, dosage_sensitive_gene_tb, 'orange', FALSE)),
                cancer_gene == 'yes' ~ map2_chr(gene_name, CNVtype, \(g, c) format_hotspots_to_badge(g,c, cancer_gene_tb, 'orange', FALSE)),
                TRUE ~ gene_name
            ),
        ) %>%
        arrange(stemcell_hotspot, cancer_gene, desc(direct_hit), start) %>%
        select(
            gene_name, gene_id, seqnames, start, end, strand, stemcell_hotspot, dosage_sensitive_gene, cancer_gene,
            any_of(extra_cols), REEV, GTex, NCBI, Ensembl
        )
    if (params$out_format == 'html') {
        colors1 <- ifelse(tb$stemcell_hotspot == 'yes', 'red' , 'white')
        colors2 <- ifelse(tb$dosage_sensitive_gene == 'yes', 'orange', 'white')
        colors3 <- ifelse(tb$cancer_gene == 'yes', 'orange', 'white')
        dt <- datatable(
            tb, rownames = FALSE, escape = FALSE,
            options = list(
                dom = 'Bftilp', pageLength = 10,
                extensions = c('Buttons'),
                buttons = c('colvis', 'copy', 'csv', 'excel', 'print'),
                columnDefs = list(list(targets = 1:8, visible = FALSE))
                )
            ) %>%
            formatStyle('stemcell_hotspot', backgroundColor = styleRow(1:nrow(tb), colors1)) %>% # textAlign = 'center'
            formatStyle('dosage_sensitive_gene', backgroundColor = styleRow(1:nrow(tb), colors2)) %>%
            formatStyle('cancer_gene', backgroundColor = styleRow(1:nrow(tb), colors3))
        return(dt)
    } else {
        tb <- tb %>% select(gene_name, gene_id, stemcell_hotspot, dosage_sensitive_gene, cancer_gene, any_of(extra_cols))
        return(kable(tb, caption = caption))
    }
}

hotspot_table_output <- function(
    hotspots, cnv_type, plotsection, stemcell_hotspot_tb, dosage_sensitive_gene_tb, cancer_gene_tb, report_config, out_format, caption = NULL
){
    tb <- bind_rows(
        stemcell_hotspot_tb,
        dosage_sensitive_gene_tb,
        cancer_gene_tb
    ) %>%
        filter(hotspot %in% hotspots & call_type %in% c('any', cnv_type))
    
    if (out_format == 'html') {
        dt <- datatable(
            tb %>%
                select(-description) %>%
                dplyr::rename(description = description_htmllinks) %>%
                select(hotspot, call_type, list_name, description, check_score, any_of(colnames(tb))) %>%
                mutate(
                    hotspot = case_when(
                        list_name == unique(stemcell_hotspot_tb$list_name) ~ map2_chr(hotspot, call_type, \(g, c) format_hotspots_to_badge(g,c, stemcell_hotspot_tb,'red', FALSE)),
                        list_name == unique(dosage_sensitive_gene_tb$list_name) ~ map2_chr(hotspot, call_type, \(g, c) format_hotspots_to_badge(g,c, dosage_sensitive_gene_tb, 'orange', FALSE)),
                        list_name == unique(cancer_gene_tb$list_name) ~ map2_chr(hotspot, call_type, \(g, c) format_hotspots_to_badge(g,c, cancer_gene_tb, 'orange', FALSE))
                    ),
                    description = str_replace_all(description, '&#013;', '<br/>')
                ) %>%
                rename_with(format_column_names),
            caption = caption,
            options = list(
                dom = 't', 
                #extensions = c('Buttons'),
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


SNV_table_output <- function (
    SNV_table, 
    roi_gr, snv_annotation_tb, gene_snv_hotspots, mutation_snv_hotspots,
    config, report_config, out_format = 'html', caption = NULL
){

    get_color <- function(crit_reason) {
        ifelse(
            crit_reason %in% report_config$SNV_analysis$critical_reasons_with_red_highlight,
            'red',
            'orange'
        )
    }
    
    tb <- SNV_table %>%
        dplyr::rename(Chromosome = seqnames, Position = start) %>%
        mutate(
            SNV = paste0(Chromosome, ':', format(Position, big.mark='.', decimal.mark = ','), ':', REF, '>', ALT),
        )

    if (out_format == 'html') {
        
        column_help_text <- c(
            'Chromosome of the SNV/SNP',
            'Position of the SNV/SNP',
            'SNV/SNP in the format "Chr:Pos:REF>ALT"',
            'ID of the SNV/SNP from the illumina array.\\nNote: rsIDs may not always be reliable',
            #FIXME: pull from label_name_definitions.yaml
            'Designation label for the SNV/SNP (critical, unrelaible critical, protein changing, reference_genotype)',
            'Reason for the critical designation label (see config:settings:SNV_analysis:critical_SNV)',
            'Reference allele of the SNV/SNP',
            'Alternative allele of the SNV/SNP',
            'Genotype of the SNV/SNP for 2 Allels: 0 stands for the reference allele, 1 for the alternative allele',
            'Reference sample genotype of the SNV/SNP for 2 Allels (0 - allele, 1 - alternative allele)',
            'Regions of interest overlapping with this SNV/SNP',
            'Gene name affected by the SNV/SNP',
            'Effect/Annotation of the SNV/SNP on the gene (from mehari)',
            'Ensembl ID of the selected transcript for the affected gene',
            'Description of the selected transcript.\\nNote: ManeSelect designates the (medically) primary transcript of a gene.',
            'HGVS.c notation of the mutation/effect of the SNV/SNP on the transcript',
            'HGVS.p notation of the mutation/effect of the SNV/SNP on the protein',
            'GenTrain Score of the SNV/SNP',
            'GenCall Score of the SNV/SNP',
            'GenCall Score of the SNV/SNP in the reference sample'            
        )
        
        dt <- tb %>%
            mutate(
                across(c(Chromosome, SNV_label, critical_reason), as.factor),
                ROI_hits = ifelse(
                    critical_reason == 'ROI-match' & !is.na(critical_reason),
                    map_chr(
                        ROI_hit,
                        \(x) format_hotspots_to_badge(x, 'any', as_tibble(roi_gr), get_color('ROI-match'))
                    ),
                    ROI_hits
                ),
                HGVS.p = ifelse(
                    critical_reason == 'hotspot-match' & !is.na(critical_reason),
                    map_chr(
                        HGVS.p,
                        \(x) format_hotspots_to_badge(x, 'any', mutation_snv_hotspots, get_color('hotspot-match'))
                    ),
                    HGVS.p
                ),
                gene_name = ifelse(
                    critical_reason == 'hotspot-gene' & !is.na(critical_reason),
                    map_chr(
                        gene_name,
                        \(x) format_hotspots_to_badge(x, 'any', gene_snv_hotspots, get_color('hotspot-gene'))
                    ),
                    gene_name
                ),
                Annotation = ifelse(
                    critical_reason == 'critical-annotation' & !is.na(critical_reason),
                    map_chr(
                        str_replace_all(Annotation, '&', '|'), 
                        \(x) format_hotspots_to_badge(x, 'any', snv_annotation_tb, get_color('critical-annotation'))
                    ),
                    Annotation
                ) %>%
                    str_replace_all('\\|', ', ')
            ) %>%
            select(
                # 0-5
                Chromosome, Position, SNV, ID, SNV_label, critical_reason, 
                # 6-9
                REF, ALT, GT, ref_GT,
                # 10-16
                ROI_hits, gene_name, Annotation, Transcript_ID, Transcript_BioType, HGVS.c, HGVS.p,
                # 17-19
                GenTrain_Score, GenCall_Score, ref_GenCall_Score
            ) %>%
            rename_with(format_column_names) %>%
            datatable(
                rownames = FALSE,
                escape = FALSE,
                extensions = c('Buttons', 'Scroller'),
                filter = 'top',
                caption = caption,
                options = list(
                    scrollY = 300, scrollCollapse = TRUE, scrollX =  TRUE, scroller = TRUE,
                    dom = 'Bftilp',
                    buttons = c('colvis', 'copy', 'csv', 'excel', 'print'),
                    #FIXME: column width restrictions don't work properly
                    # autoWidth = TRUE, # necessary for columnDefs/width to work
                    # fixedColumns = TRUE,
                    columnDefs = list(
                        #This uses 0-indexing vs the usual R 1-indexing
                        list(targets = c(0:1,3,5:7,13:15,17,19), visible = FALSE)
                        # list(targets = c(12), width = '300px'),
                        # list(targets = c(15,16), width = '200px')
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
        return(dt)
    } else {
        tb <- tb %>%
            select(
                SNV, SNV_label, critical_reason, GT, ref_GT,
                ROI_hits, gene_name, Annotation, HGVS.p,
                GenCall_Score
            ) %>%
            rename_with(format_column_names)
        return(kable(tb, caption = caption))
    }
}