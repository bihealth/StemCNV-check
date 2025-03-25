suppressMessages(library(tidyverse))

minimal_hotspots <- tibble::tibble(
    list_name = 'test-list',
    hotspot = c('DDX11L1', 'dummyC', 'chr1:40000-50000', '1p36', '1p35.2', '1q21'),
    mapping = c('gene_name', 'gene_name', 'position', 'gband', 'gband', 'gband'),
    call_type = c('any', 'gain', 'any', 'loss', 'any', 'gain'),
    check_score = c(30, 15, 10, 10, 12, 10),
    description = c('Sources: dummy', 'Something: Dummy{1}', 'dummy', 'Sources: dummy{1}\\nSomething: else{2}',
                    'dummy', 'Sources: dummy{1},dummy{2}'),
    description_doi = c(NA, 'doi', NA, 'doi, doi2', NA, 'doi, doi2'),
    description_htmllinks = c(
        'Sources: dummy', 
        'Something: <a href="https://doi.org/doi" target="_blank" rel="noopener noreferrer">Dummy</a>',
        'dummy',
        'Sources: <a href="https://doi.org/doi" target="_blank" rel="noopener noreferrer">dummy</a>&#013;Something: <a href="https://doi.org/doi2" target="_blank" rel="noopener noreferrer">else</a>',
        'dummy', 
        'Sources: <a href="https://doi.org/doi" target="_blank" rel="noopener noreferrer">dummy</a>,<a href="https://doi.org/doi2" target="_blank" rel="noopener noreferrer">dummy</a>'
    )
)

minimal_hotspots_positions <- minimal_hotspots %>%
    mutate(
        seqnames = 'chr1',
        start = c(11873, 28050000, 40000, 0, 30200000, 142600000),
        end = c(14409, 28070000, 50000, 7200000, 32400000, 155000000),
        strand = c('+', '+', '*', '*', '*', '*'),
    )

meahri_ann_header <- c(
    'Allele', 'Annotation', 'Annotation_Impact', 'gene_name', 'gene_id', 'Feature_Type', "Feature_ID",
    'Transcript_BioType', 'Rank', 'HGVS.c', 'HGVS.p', 'cDNA.pos', 'CDS.pos', 'AA.pos',        
    # cDNA.pos / cDNA.length | CDS.pos / CDS.length | AA.pos / AA.length |
    'Distance', 'Strand', 'errors' # ERRORS / WARNINGS / INFO
)
# Directly setup the expected SNV table from vcf parsing
sample_SNV_tb <- suppressMessages(read_tsv(test_path('../data/minimal_probes_annotated.vcf'), skip = 8)) %>%
    dplyr::rename(
        seqnames = `#CHROM`,
        start = POS,
        sample_id = FORMAT
    ) %>%
    mutate(
        end = start,
        width = 1,
        strand = '*',
        sample_id = 'annot_sample',
        INFO = ifelse(INFO == '.', ';', INFO) %>%
            str_replace('^ANN=', ';ANN='),
        ROI_hits = NA_character_
    ) %>%
    separate(INFO, into = c('GenTrain_Score', 'ANN'), sep = ';', fill='right') %>%
    separate(annot_sample, into = c('GT', 'IGC', 'LRR', 'BAF'), sep = ':') %>%
    mutate(
        ANN = str_replace(ANN, '^ANN=', ''),
        GenTrain_Score = str_replace(GenTrain_Score, '^GenTrain_Score=', '') %>% as.numeric(),
        GT = ifelse(GT == '.', './.', GT),
        across(c(ID, QUAL), ~ifelse(. == '.', NA_character_, .))
    ) %>%
    separate(ANN, into = meahri_ann_header, sep = '\\|', fill = 'right') %>%
    relocate(seqnames, start, end, width, strand) %>%
    dplyr::rename(GenCall_Score = IGC, Transcript_ID  = Feature_ID) %>%
    select(-LRR, -BAF) %>%
    mutate(
        Allele = ifelse(Allele == '', NA_character_, Allele),
        across(everything(), ~as.character(.)),
        across(c(start, end, width, QUAL), ~as.numeric(.)),
        seqnames = factor(seqnames),
        strand = factor(strand, levels = c('+', '-', '*')),
    )
