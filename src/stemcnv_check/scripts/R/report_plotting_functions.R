library(tidyverse)
library(patchwork)
library(ggrepel)

add_CNV_plot_styling <- function (gg, panel_space_val = unit(5, units = 'mm')) {
    gg +
        facet_wrap(~Sample_ID, nrow = 1) +
        theme(
            strip.background = element_blank(),
            strip.text.x = element_blank(),
            panel.spacing = panel_space_val
        )
}

make_BAF_panel <- function(
    chr, win_start, win_end, 
    plot.data,
    highlight.tb = NULL, 
    area_tb = NULL, area_alpha = 0.3
) {    
    gg <- plot.data %>%
        mutate(color = ifelse(filter.passed, 'blue', 'grey70')) %>%
        ggplot() +
        geom_hline(yintercept = 0, col = 'black', linewidth=0.5) +
        geom_hline(yintercept = 1, col = 'black', linewidth=0.5)
    if(!is.null(area_tb)){
        gg <- gg + geom_rect(
            data = area_tb,
            aes(xmin = start, xmax = end, fill = color,ymin = 0, ymax = 1),
            alpha = area_alpha
        ) +
            scale_fill_identity()
    }
    gg <- gg +
        geom_point(
            aes(x = Position, y = `B Allele Freq`, color = color),
            size = 0.5, shape = 20, show.legend = F
        ) 
    if (!is.null(highlight.tb)) {
        highlight.plotdata <- left_join(
            highlight.tb %>% select(Sample_ID, ID, color),
            plot.data %>% select(Sample_ID, ID, Position, `B Allele Freq`),
            by = c('Sample_ID', 'ID')
        )        
        gg <- gg + geom_point(
            data = highlight.plotdata,
            aes(x = Position, y = `B Allele Freq`, color = color),
            size = 1, shape = 20, show.legend = F
        )
    }    
    gg +
        scale_color_identity() +
        theme_classic() +
        scale_x_continuous(
            expand = expansion(),
            labels = label_number(big.mark = '', decimal.mark = ','),
            limits = c(win_start, win_end),
            oob = oob_keep
        ) +
        scale_y_continuous(expand = expansion(), limits = c(-0.1, 1.1), oob = oob_squish, breaks = c(0, 0.5, 1)) +
        labs(y = 'B Allele Frequency', x = paste0('Position (', chr, ')')) 
}

make_LRR_panel <- function(
    chr, win_start, win_end, 
    plot.data,
    highlight.tb = NULL, 
    area_tb = NULL, area_alpha = 0.3
) {    
    gg <- plot.data %>%
        mutate(color = ifelse(filter.passed, 'blue', 'grey70')) %>%
        ggplot()
    if(!is.null(area_tb)){
        gg <- gg + geom_rect(
            data = area_tb,
            aes(xmin = start, xmax = end, fill = color,ymin = -1.5, ymax = 1.5),
            alpha = area_alpha
        ) +
            scale_fill_identity()
    }
    gg <- gg +
        geom_hline(yintercept = 0, col = 'grey10', linewidth=0.2) +
        geom_point(
            aes(x = Position, y = `Log R Ratio`,color = color),
            size = 0.5, shape = 20, show.legend = F
        ) 
    if (!is.null(highlight.tb)) {
        highlight.plotdata <- left_join(
            highlight.tb %>% select(Sample_ID, ID, color),
            plot.data %>% select(Sample_ID, ID, Position, `Log R Ratio`),
            by = c('Sample_ID', 'ID')
        )        
        gg <- gg + geom_point(
            data = highlight.plotdata,
            aes(x = Position, y = `Log R Ratio`, color = color),
            size = 1, shape = 20, show.legend = F
        )
    }    
    gg +
        scale_color_identity() +
        theme_classic() +
        scale_x_continuous(
            expand = expansion(),
            labels = label_number(big.mark = '', decimal.mark = ','),
            limits = c(win_start, win_end),
            position = 'top'
        ) +
        scale_y_continuous(expand = expansion(), limits = c(-1.5, 1.5), oob = oob_squish) +
        labs(y = 'Log R Ratio', x = paste0('Position (', chr, ')')) +
        facet_wrap(~Sample_ID, nrow = 1)
}

make_CNV_panel <- function(
    chr, win_start, win_end, 
    calls,
    label_column = 'call_label'
) {
    ggplot(calls) +
        geom_tile(aes(x = x_pos, y = y_pos, width = width, height = .9, fill = color)) +
        scale_fill_identity() +
        geom_text(
            aes(label = !!sym(label_column), x = x_pos, y = y_pos),
            vjust = 0.5, hjust = 0.5, size = 2.5
        ) +
        scale_x_continuous(
            expand = expansion(),
            labels = label_number(big.mark = '', decimal.mark = ','),
            limits = c(win_start, win_end),
            oob = oob_keep
        ) +
        scale_y_continuous(expand = expansion()) +
        theme_classic() +
        theme(
            axis.title = element_blank(),
            axis.line = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            plot.background = element_blank(),
            panel.background = element_blank(),
            panel.border = element_blank()
        ) +
        labs(y = 'Calls') # x = paste0('Position (', chr, ')')
}

make_gene_data <- function(
    chr, win_start, win_end, sample_headers,
    gr_genes, area_tb,
    stemcell_hotspot_list, dosage_sensitive_gene_list, cancer_gene_list    
) {
    direct_genes <- gr_genes %>%
        filter_by_overlaps(GRanges(seqnames = chr, strand = '*', ranges = IRanges(start = area_tb$start, end = area_tb$end))) %>%
        as_tibble()
    
    gr_genes %>%
        filter_by_overlaps(GRanges(seqnames = chr, strand = '*', ranges = IRanges(start = win_start, end = win_end))) %>%
        as_tibble() %>%
        mutate(
            x_pos = (end + start) / 2,
            y_pos = ifelse(strand == '+', 1, 0),
            Sample_ID = paste(sample_headers, collapse = '---'),
            direct_hit = gene_id %in% direct_genes$gene_id,
            stemcell_hotspot = gene_name %in% stemcell_hotspot_list,
            dosage_sensitive_gene = gene_name %in% dosage_sensitive_gene_list,
            cancer_gene = gene_name %in% cancer_gene_list,
        ) %>%
        separate_rows(Sample_ID, sep = '---') %>%
        # Need to ensure table contains all samples (reference) so everything is properly facet_wrapped
        bind_rows(
            tibble(
                Sample_ID = factor(sample_headers, levels = sample_headers),
                x_pos = NA_integer_, y_pos = NA_integer_
            )
        )
    
}

make_gene_panel <- function(
    chr, win_start, win_end,
    gene.data
) {
    
    gene_track <- ggplot(gene.data) +
        geom_tile(
            aes(
                x = x_pos, y = y_pos, width = width, height = .9,
                fill = case_when(
                    stemcell_hotspot ~ 'red',
                    dosage_sensitive_gene ~ 'orange',
                    cancer_gene   ~ 'orange',
                    direct_hit  ~ 'black',
                    TRUE        ~ 'grey50'
                )
            ),
            show.legend = F
        ) +
        scale_x_continuous(expand = expansion(), limits = c(win_start, win_end), oob = oob_keep) +
        scale_y_continuous(expand = expansion(add = c(0.25, 0.25))) +
        scale_fill_identity() +
        theme_void() +
        theme(
            axis.title.y = element_text(angle = 90, vjust = 1)
        ) +
        labs(y = 'Genes') # x = paste0('Position (', chr, ')')
    
}

make_header_data <- function(
    chr, win_start, win_end, sample_headers,
    gr_info, stemcell_hotspot_list
) {
   
    info_data <- gr_info %>%
        filter_by_overlaps(GRanges(seqnames = chr, strand = '*', ranges = IRanges(start = win_start, end = win_end))) %>%
        as_tibble() %>%
        mutate(
            x_pos = (end + start) / 2,
            y_pos = 0,
            Sample_ID = paste(sample_headers, collapse = '---'),
            color = case_when(
                str_detect(section_name, paste(stemcell_hotspot_list, collapse = '|')) ~ 'red',
                # Onlt the stemcell list contains gbands
                # str_detect(section_name, paste(cancer_gene_list, collapse = '|'))   ~ 'orange',
                band_staining == 'gpos100' ~ 'black',
                band_staining ==  'gpos50' ~ 'grey30',
                band_staining ==  'gpos25' ~ 'grey70',
                centromer ~ 'lightblue',
                TRUE ~ 'white'
            ),
            textcolor = ifelse(color %in% c('black', 'grey30'), 'white', 'black')
        ) %>%
        separate_rows(Sample_ID, sep = '---') %>%
        mutate(Sample_ID = factor(Sample_ID, levels = sample_headers)) %>%
        # Need to ensure table contains reference so everything is properly facet_wrapped
        bind_rows(
            tibble(
                Sample_ID = factor(sample_headers, levels = sample_headers),
                x_pos = NA_integer_, y_pos = NA_integer_
            )
        )  # +
        # labs(y = 'gBand', x = paste0('Position (', chr, ')'))
}

make_header_panel <- function(
    chr, win_start, win_end, header_data, label_colunm = 'section_name'
) {
    
    ggplot(header_data) + 
        theme_classic() +
        geom_tile(
            aes(x = x_pos, y = y_pos, width = width, height = .9, fill = color),
            color = 'black', linewidth = 0.2
        ) +
        scale_fill_identity() +
        # Use repel to keep gband names in the plot area
        geom_text(
            aes(label = !!sym(label_colunm), x = x_pos, y = y_pos, color = textcolor),
            vjust = 0.5, hjust = 0.5, size = 2.5, show.legend = F
        ) +
        scale_color_identity() +
        theme_classic() +
        scale_x_continuous(
            expand = expansion(),
            labels = label_number(big.mark = '', decimal.mark = ','),
            limits = c(win_start, win_end),
            oob = oob_keep
        ) +
        scale_y_continuous(expand = expansion()) +
        theme(
            axis.line = element_blank(),
            axis.title = element_blank(),
            axis.ticks = element_blank(),
            axis.text = element_blank(),
            plot.background = element_blank(),
            panel.background = element_blank(),
            panel.border = element_blank()
        )
}


make_call_plot <- function(
    call.row, raw_LRR_BAF, cnv_calls, gr_genes, gr_info,
    sample_headers,
    total_min_size = 2e6, flank_factor = 2
) {
    
    chr <- call.row$chrom
    area_tb <- tibble(
        Sample_ID = factor(sample_headers, levels = sample_headers),
        start = call.row$start, end = call.row$end, color = 'grey70'        
    )
    
    if (total_min_size < 1e4) {
        warning(str_glue('Re-Setting Minimum size window around primary plot area (CNV/region: {chr}:{area_tb$start}-{area_tb$end}) to at least 10kb'))
        total_min_size <- 1e4
    }
    
    # Set initial plot window as Call size * (1+2*flank_factor)
    win_start <- call.row$start - call.row$Size * flank_factor
    win_end <- call.row$end + call.row$Size * flank_factor
    # If it is below minimum increase accordingly
    if (win_end - win_start < total_min_size) {
        #Don't want .5 positions
        win_start <- win_start - floor((total_min_size - win_end + win_start)/2)
        win_end <- win_end + ceiling((total_min_size - win_end + win_start)/2)
    }
    
    #Shrink down ends if they overlap 0/Chr_Max (from what is covered by probes)
    chr_max <- raw_LRR_BAF %>% filter(Chr == chr) %>% pull(Position) %>% max()
    win_start <- max(win_start, 0)
    win_end <- min(win_end, chr_max)

    # get raw LRR & BAF data; mark filtered points
    plot.data <- raw_LRR_BAF %>%
        # Assume all samples in sample_headers have been loaded
        filter(sample_id %in% names(sample_headers) & Chr == chr & Position >= win_start & Position <= win_end) %>%
        mutate(Sample_ID = factor(sapply(sample_id, function(x) sample_headers[x]), levels = sample_headers))
    
    if (nrow(plot.data)==0){
        warn_msg <- str_glue('No SNP probes found in primary plot area: {chr}:{win_start}-{win_end}')
        warning(warn_msg)
        return(list('gg' = warn_msg, 'genes' = tibble(), 'hotspots' = c()))
    }

    cnv_tools <- cnv_calls$CNV_caller %>% unique() %>% str_subset('StemCNV-check', TRUE) %>% sort()
    get_cnv_y <- function(CNV_type, CNV_caller) {
        # Go by tool order (start from 1)
        out <- match(CNV_caller, cnv_tools)
        # For PennCNV LOH & gain/loss may overlap and need separate tracks
        out <- ifelse(CNV_type %!in% c('gain', 'loss'), 1-out, out)
        out
    }
    
    calls <- cnv_calls %>%
        filter(sample_id %in% names(sample_headers) & seqnames == chr & end >= win_start & start < win_end) %>%
        as_granges() %>%
        unsplit_merged_CNV_callers() %>%
        as_tibble() %>% 
        mutate(
            x_pos = (end + start) / 2,
            y_pos = map2_int(CNV_type, CNV_caller, get_cnv_y),
            color = case_when(
                CNV_type == 'gain' ~ '#1a9850',
                CNV_type == 'loss' ~ '#f46d43',
                CNV_type == 'LOH'  ~ 'grey50'
            ),
            Sample_ID = factor(sapply(sample_id, function(x) sample_headers[x]), levels = sample_headers),
            call_label = str_glue('{CNV_caller}: {CNV_type}')
        ) %>%
        # Need to ensure table contains reference so everything is properly facet_wrapped
        bind_rows(
            tibble(
                Sample_ID = factor(sample_headers, levels = sample_headers),
                x_pos = NA_integer_, y_pos = NA_integer_
            )
        )
    
    stemcell_hotspot_list <- call.row$stemcell_hotspot %>% str_split('\\|') %>% unlist()
    dosage_sensitive_gene_list <- call.row$dosage_sensitive_gene %>% str_split('\\|') %>% unlist()
    cancer_gene_list <- call.row$cancer_gene %>% str_split('\\|') %>% unlist()
    panel_space_val <- unit(5, units = 'mm')
    header.data <- make_header_data(
        chr, win_start, win_end, sample_headers, 
        gr_info, stemcell_hotspot_list
    )
    gene.data <- make_gene_data(
        chr, win_start, win_end, sample_headers,
        gr_genes, area_tb,
        stemcell_hotspot_list, dosage_sensitive_gene_list, cancer_gene_list 
    )
        
    header <- make_header_panel(chr, win_start, win_end, header.data) +
        facet_wrap(~Sample_ID, nrow = 1)
    cnv_track <- make_CNV_panel(chr, win_start, win_end, calls) %>% 
        add_CNV_plot_styling(panel_space_val)
    lrr <- make_LRR_panel(
        chr, win_start, win_end, plot.data,
        area_tb = area_tb 
    ) %>% 
        add_CNV_plot_styling(panel_space_val)
    baf <- make_BAF_panel(
        chr, win_start, win_end, plot.data,
        area_tb = area_tb
    ) %>%
        add_CNV_plot_styling(panel_space_val)
    gene_track <- make_gene_panel(chr, win_start, win_end, gene.data) %>%
        add_CNV_plot_styling(panel_space_val)

    n_cnvs <- length(na.omit(unique(calls$CNV_type)))
    gg <- header / cnv_track / lrr / baf / cnv_track / gene_track + plot_layout(heights = c(1, n_cnvs, 10, 10, n_cnvs, 2))

    gene.data <- gene.data %>%
        filter(!is.na(x_pos)) %>%
        dplyr::select(
            seqnames, start, end, width, strand, stemcell_hotspot, dosage_sensitive_gene, cancer_gene, 
            direct_hit, gene_name, gene_type, gene_id
        ) %>%
        mutate(CNVtype = as.character(call.row$CNV_type)) %>%
        unique()

    list(
        'gg' = gg,
        'genes' = gene.data,
        'hotspots' = c(stemcell_hotspot_list, dosage_sensitive_gene_list, cancer_gene_list) %>% na.omit()
    )

}

make_chromsome_overview_plot <- function(
    chr, sample_headers, label_SNVs,
    raw_LRR_BAF, cnv_calls, gr_info
) {
    plot.data <- raw_LRR_BAF %>%
        filter(sample_id %in% names(sample_headers) & Chr == chr ) %>%
        mutate(Sample_ID = factor(sapply(sample_id, function(x) sample_headers[x]), levels = sample_headers))
    
    area_tb <- cnv_calls %>%
        filter(sample_id %in% names(sample_headers) & seqnames == chr) %>%
        mutate(
            color = case_when(
                CNV_type == 'gain' ~ '#1a9850',
                CNV_type == 'loss' ~ '#f46d43',
                CNV_type == 'LOH'  ~ 'grey50'
            ),
            Sample_ID = factor(sapply(sample_id, function(x) sample_headers[x]), levels = sample_headers)
        ) %>%
        # Need to ensure table contains reference so everything is properly facet_wrapped
        bind_rows(
            tibble(
                Sample_ID = factor(sample_headers, levels = sample_headers)
            )
        )
    
    SNV_label_data <- label_SNVs %>%
        filter(sample_id %in% names(sample_headers) & seqnames == chr) %>%
        mutate(
            Sample_ID = factor(sapply(sample_id, function(x) sample_headers[x]), levels = sample_headers),
        ) %>%
        # Need to ensure table contains reference so everything is properly facet_wrapped
        bind_rows(
            tibble(
                Sample_ID = factor(sample_headers, levels = sample_headers),
            )
        )
    
    chr_end <- gr_info %>% as_tibble() %>% filter(seqnames == chr) %>% pull(end) %>% max()
    
    header.data <- make_header_data(
        #need the full data here
        chr = chr, win_start = 0, win_end = chr_end, sample_headers = sample_headers, 
        gr_info = gr_info, stemcell_hotspot_list = 'dummy_string'
    ) %>%
        mutate(dummy_label = NA_character_)
    panel_space_val <- unit(5, units = 'mm')
    
    header <- make_header_panel(chr, 0, chr_end, header.data, 'dummy_label')  +
        facet_wrap(~Sample_ID, nrow = 1)
    lrr <- make_LRR_panel(chr, 0, chr_end, plot.data, SNV_label_data, area_tb, 0.7) %>% 
        add_CNV_plot_styling(panel_space_val)
    baf <- make_BAF_panel(chr, 0, chr_end, plot.data, SNV_label_data, area_tb, 0.7) %>%
        add_CNV_plot_styling(panel_space_val)
    
    gg <- header / lrr / baf  + plot_layout(heights = c(1, 8, 8))
    
    gg
}