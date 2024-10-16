library(tidyverse)
library(patchwork)
library(ggrepel)


make_LRR_BAF_plots <- function(
    call.row, raw_LRR_BAF, cnv_calls, gr_genes, gr_info,
    total_min_size = 2e6, flank_factor = 2
) {
    
    chr <- call.row$chrom
    
    if (total_min_size < 1e4) {
        warning(str_glue('Re-Setting Minimum size window around primary plot area (CNV/region: {chr}:{call.row$start}-{call.row$end}) to at least 10kb'))
        total_min_size <- 1e4
    }
    
    cnv_tools <- cnv_calls$CNV_caller %>% unique() %>% str_subset('StemCNV-check', TRUE) %>% sort()
    get_cnv_y <- function(CNV_type, CNV_caller) {
        # Go by tool order (start from 1)
        out <- match(CNV_caller, cnv_tools)
        # For PennCNV LOH & gain/loss may overlap and need separate tracks
        out <- ifelse(CNV_type %!in% c('gain', 'loss'), 1-out, out)
        out
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
        filter(Chr == chr & Position >= win_start & Position <= win_end) %>%
        mutate(Sample_Name = factor(sapply(sample_id, function(x) sample_headers[x]), levels = sample_headers))
    
    if (nrow(plot.data)==0){
        warn_msg <- str_glue('No SNP probes found in primary plot area: {chr}:{win_start}-{win_end}')
        warning(warn_msg)
        return(list('gg' = warn_msg, 'genes' = tibble(), 'hotspots' = c()))
    }

    calls <- cnv_calls %>%
        filter(seqnames == chr & end >= win_start & start < win_end) %>%
        as_granges() %>%
        unsplit_merged_CNV_callers() %>%
        as_tibble() %>% 
        mutate(
            x_pos = (end + start) / 2,
            y_pos = map2_int(CNV_type, CNV_caller, get_cnv_y),
            color = ifelse(CNV_type %!in% c('gain', 'loss'), 'grey50', '#1a9850'),
            color = ifelse(CNV_type == 'loss', '#f46d43', color),
            Sample_Name = factor(sapply(sample_id, function(x) sample_headers[x]), levels = sample_headers)
        ) %>%
        # Need to ensure table contains reference so everything is properly facet_wrapped
        bind_rows(
            tibble(
                Sample_Name = factor(sample_headers, levels = sample_headers),
                x_pos = NA_integer_, y_pos = NA_integer_
            )
        )

    direct_genes <- gr_genes %>%
        filter_by_overlaps(GRanges(seqnames = chr, strand = '*', ranges = IRanges(start = call.row$start, end = call.row$end))) %>%
        as_tibble()

    high_impact_list <- call.row$HighImpact %>% str_split('\\|') %>% unlist()
    highlight_list <- call.row$Highlight %>% str_split('\\|') %>% unlist()
    gene.data <- gr_genes %>%
        filter_by_overlaps(GRanges(seqnames = chr, strand = '*', ranges = IRanges(start = win_start, end = win_end))) %>%
        as_tibble() %>%
        mutate(
            x_pos = (end + start) / 2,
            y_pos = ifelse(strand == '+', 1, 0),
            Sample_Name = paste(sample_headers, collapse = '---'),
            direct_hit = gene_id %in% direct_genes$gene_id,
            high_impact = gene_name %in% high_impact_list,
            highlight = gene_name %in% highlight_list,
        ) %>%
        separate_rows(Sample_Name, sep = '---') %>%
        # Need to ensure table contains reference so everything is properly facet_wrapped
        bind_rows(
            tibble(
                Sample_Name = factor(sample_headers, levels = sample_headers),
                x_pos = NA_integer_, y_pos = NA_integer_
            )
        )

    info_data <- gr_info %>%
        filter_by_overlaps(GRanges(seqnames = chr, strand = '*', ranges = IRanges(start = win_start, end = win_end))) %>%
        as_tibble() %>%
        mutate(
            x_pos = (end + start) / 2,
            y_pos = 0,
            Sample_Name = paste(sample_headers, collapse = '---'),
            color = case_when(
                str_detect(section_name, paste(high_impact_list, collapse = '|')) ~ 'red',
                str_detect(section_name, paste(highlight_list, collapse = '|'))   ~ 'orange',
                band_staining == 'gpos100' ~ 'black',
                band_staining ==  'gpos50' ~ 'grey30',
                band_staining ==  'gpos25' ~ 'grey70',
                centromer ~ 'lightblue',
                TRUE ~ 'white'
            ),
            textcolor = ifelse(color %in% c('black', 'grey30'), 'white', 'black')
        ) %>%
        separate_rows(Sample_Name, sep = '---') %>%
        mutate(Sample_Name = factor(Sample_Name, levels = sample_headers)) %>%
        # Need to ensure table contains reference so everything is properly facet_wrapped
        bind_rows(
            tibble(
                Sample_Name = factor(sample_headers, levels = sample_headers),
                x_pos = NA_integer_, y_pos = NA_integer_
            )
        )

    panel_space_val <- unit(5, units = 'mm')

    cnv_track <- ggplot(calls) +
        geom_tile(aes(x = x_pos, y = y_pos, width = width, height = .9, fill = color)) +
        scale_fill_identity() +
        geom_text(
            aes(label = paste0(CNV_caller, ': ', CNV_type), x = x_pos, y = y_pos),
            vjust = 0.5, hjust = 0.5, size = 2.5
        ) +
        scale_x_continuous(
            expand = expansion(),
            labels = label_number(big.mark = '.', decimal.mark = ','),
            limits = c(win_start, win_end),
            oob = oob_keep
        ) +
        scale_y_continuous(expand = expansion()) +
        facet_wrap(~Sample_Name, nrow = 1) +
        theme_classic() +
        theme(
            axis.title = element_blank(),
            axis.line = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            plot.background = element_blank(),
            panel.background = element_blank(),
            panel.border = element_blank(),
            strip.background = element_blank(),
            strip.text.x = element_blank(),
            panel.spacing = panel_space_val
        ) +
        labs(y = 'Calls')

    gene_track <- ggplot(gene.data) +
        geom_tile(
            aes(
                x = x_pos, y = y_pos, width = width, height = .9,
                fill = case_when(
                    high_impact ~ 'red',
                    highlight   ~ 'orange',
                    direct_hit  ~ 'black',
                    TRUE        ~ 'grey50'
                )
            ),
            show.legend = F
        ) +
        scale_x_continuous(expand = expansion(), limits = c(win_start, win_end), oob = oob_keep) +
        scale_y_continuous(expand = expansion(add = c(0.25, 0.25))) +
        scale_fill_identity() +
        facet_wrap(~Sample_Name, nrow = 1) +
        theme_void() +
        theme(
            strip.background = element_blank(),
            strip.text.x = element_blank(),
            axis.title.y = element_text(angle = 90, vjust = 1),
            panel.spacing = panel_space_val
        ) +
        labs(y = 'Genes')

    lrr <- ggplot(plot.data) +
        geom_rect(
            data = tibble(Sample_Name = sample_headers),
            aes(xmin = call.row$start, xmax = call.row$end, ymin = -1.5, ymax = 1.5),
            fill = 'grey50', alpha = 0.3
        ) +
        geom_hline(yintercept = 0, col = 'grey10', linewidth=0.2) +
        geom_point(
            aes(x = Position, y = `Log R Ratio`,color = filter.passed),
            size = 0.5, shape = 20, show.legend = F
        ) +
        scale_color_manual(values=c('TRUE' = 'blue', 'FALSE' = 'grey70')) +
        theme_classic() +
        scale_x_continuous(
            expand = expansion(),
            labels = label_number(big.mark = '.', decimal.mark = ','),
            limits = c(win_start, win_end),
            position = 'top'
        ) +
        scale_y_continuous(expand = expansion(), limits = c(-1.5, 1.5), oob = oob_squish) +
        labs(y = 'Log R Ratio', x = paste0('Position (', chr, ')')) +
        facet_wrap(~Sample_Name, nrow = 1) +
        theme(
            strip.background = element_blank(),
            strip.text.x = element_blank(),
            panel.spacing = panel_space_val
        )

    baf <- ggplot(plot.data) +
        geom_hline(yintercept = 0, col = 'black', linewidth=0.5) +
        geom_hline(yintercept = 1, col = 'black', linewidth=0.5) +
        geom_rect(
            data = tibble(Sample_Name = sample_headers),
            aes(xmin = call.row$start, xmax = call.row$end, ymin = 0, ymax = 1),
            fill = 'grey50', alpha = 0.3
        ) +
        geom_point(
            aes(x = Position, y = `B Allele Freq`,color = filter.passed),
            size = 0.5, shape = 20, show.legend = F
        ) +
        scale_color_manual(values=c('TRUE' = 'blue', 'FALSE' = 'grey70')) +
        theme_classic() +
        scale_x_continuous(
            expand = expansion(),
            labels = label_number(big.mark = '.', decimal.mark = ','),
            limits = c(win_start, win_end),
            oob = oob_keep
        ) +
        scale_y_continuous(expand = expansion(), limits = c(-0.1, 1.1), oob = oob_squish, breaks = c(0, 0.5, 1)) +
        labs(y = 'B Allele Frequency', x = paste0('Position (', chr, ')')) +
        facet_wrap(~Sample_Name, nrow = 1) +
        theme(
            strip.background = element_blank(),
            strip.text.x = element_blank(),
            panel.spacing = panel_space_val
        )

    header <- ggplot(info_data) + 
        facet_wrap(~Sample_Name, nrow=1) + theme_classic() +
        geom_tile(
            aes(x = x_pos, y = y_pos, width = width, height = .9, fill = color),
            color = 'black', linewidth = 0.2
        ) +
        scale_fill_identity() +
        # Use repel to keep gband names in the plot area
        geom_text(
            aes(label = section_name, x = x_pos, y = y_pos, color = textcolor),
            vjust = 0.5, hjust = 0.5, size = 2.5, show.legend = F
        ) +
        scale_color_identity() +
        theme_classic() +
        scale_x_continuous(
            expand = expansion(),
            labels = label_number(big.mark = '.', decimal.mark = ','),
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
            panel.border = element_blank(),
            panel.spacing = panel_space_val
        )

    n_cnvs <- length(na.omit(unique(calls$CNV_type)))
    gg <- header / cnv_track / lrr / baf / cnv_track / gene_track + plot_layout(heights = c(1, n_cnvs, 10, 10, n_cnvs, 2))

    gene.data <- gene.data %>%
        filter(!is.na(x_pos)) %>%
        dplyr::select(seqnames, start, end, width, strand, high_impact, highlight, direct_hit, gene_name, gene_type, gene_id) %>%
        mutate(CNVtype = as.character(call.row$CNV_type)) %>%
        unique()

    list('gg' = gg, 'genes' = gene.data, 'hotspots' = c(high_impact_list, highlight_list) %>% na.omit())

}
