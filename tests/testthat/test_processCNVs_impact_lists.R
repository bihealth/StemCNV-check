library(testthat)

library(tidyverse)
library(plyranges)

source(test_path('../../stemcnv_check/scripts/R/helper_functions.R'))
source(test_path('../../stemcnv_check/scripts/R/processCNVs_annotate_impact_lists.R'))


config <- list(
    'snakedir' = '',
    'static_data' = list(
        'genome_gtf_file' = test_path('../data/hg_minimal.gtf'),
        'genomeInfo_file' = test_path('../data/gr_info_minimal.tsv')
    ),
    'settings' = list(
        'CNV_processing' = list(
            'gene_overlap' = list(
                'exclude_gene_type_regex' = c(),
                'include_only_these_gene_types' = c('lncRNA', 'miRNA', 'protein_coding'),
                'high_impact_list' = test_path('../data/minimal-hotspots.tsv')
            )
        )
    )
)

# Functions to test
test_that('tb_to_gr_by_position', {
    tb_pos <- tibble(
        pos_col = c('chr1:100-20000', 'chrX:30000-40000', '4:300-500'),
        metacol1 = c(1, 2, 3),
        metacol2 = list(c('A', 'B'), c(1, 2), c(5))
    )
    gr_pos <- tibble(
        seqnames = c('chr1', 'chrX', '4'),
        start = c(100, 30000, 300),
        end = c(20000, 40000, 500),
        pos_col = c('chr1:100-20000', 'chrX:30000-40000', '4:300-500'),
        metacol1 = c(1, 2, 3),
        metacol2 = list(c('A', 'B'), c(1, 2), c(5))
    ) %>% as_granges()

    expect_equal(tb_to_gr_by_position(tb_pos, 'pos_col'), gr_pos)

    tb_pos <- bind_rows(
        tb_pos,
        tibble(pos_col = 'XXX:1234')
    )
    expect_error(tb_to_gr_by_position(tb_pos, 'pos_col'))

    #test empty input; expect_qual is tricky due to info in seqnames
    expect_length(
        tb_to_gr_by_position(tb_pos %>% filter(metacol1 == 5), 'pos_col'),
        0
    )
})

#TODO: somehow this didn't/doesn't catch all possible issues
# function failed before, due to not correctly checking for all matched gband names
test_that('tb_to_gr_by_gband', {
    gr_info <- load_genomeInfo(config)
    
    tb_band <- tibble(
        band_col = c('1p35', '1p11', '1q21.2'),
        metacol1 = c(1, 2, 3),
        metacol2 = list(c('A', 'B'), c(1, 2), c(5))
    )
    gr_band <- tibble(
        seqnames = c('chr1', 'chr1', 'chr1'),
        start = c(28000000, 120600000, 147000000),
        end = c(34600000, 125000000, 150300000),
        band_col = c('1p35', '1p11', '1q21.2'),
        metacol1 = c(1, 2, 3),
        metacol2 = list(c('A', 'B'), c(1, 2), c(5))
    ) %>% as_granges()

    expect_equal(tb_to_gr_by_gband(tb_band,  gr_info, 'band_col'), gr_band)

    # Test error on wrong formatted gband
    tb_band2 <- bind_rows(
        tb_band,
        tibble(band_col = 'XXX')
    )
    expect_error(tb_to_gr_by_gband(tb_band2,  gr_info, 'band_col'))

    # Test error on non-annotated gband
    tb_band2 <- bind_rows(
        tb_band,
        tibble(band_col = '1q44')
    )
    expect_error(tb_to_gr_by_gband(tb_band2,  gr_info, 'band_col'))

    #test empty input; expect_qual is tricky due to info in seqnames
    expect_length(
        tb_to_gr_by_gband(tb_band %>% filter(metacol1 == 5), gr_info, 'band_col'),
        0
    )

})

test_that('parse_hotspot_table', {
    tb <- load_hotspot_table(config)
    gr_info <- load_genomeInfo(config)
    gr_genes <- load_gtf_data(config)

    expected_gr <- minimal_probes %>%
        mutate(
            seqnames = 'chr1',
            start = c(11873, 28050000, 40000, 0, 30200000, 142600000),
            end = c(14409, 28070000, 50000, 7200000, 32400000, 155000000),
            strand = c('+', '+', '*', '*', '*', '*'),
        ) %>%
        as_granges()

    expect_equal(parse_hotspot_table(tb, gr_genes, gr_info), expected_gr)
    
    tb <- bind_rows(tb, tibble(mapping = 'unclear'))
    expect_error(parse_hotspot_table(tb, gr_genes, gr_info))
    
    #test empty input
    expect_length(
        parse_hotspot_table(tb %>% filter(description_doi == '123'), gr_genes, gr_info),
        0
    )
})

# FIXME: enable skipping on all but manual execution
# test_that('parse inbuilt tables', {
#     config <- list(
#         'snakedir' = test_path('../../stemcnv_check/'),
#         'static_data' = list(
#             'genome_gtf_file' = test_path('../../test_folders/static-data/gencode.v42.basic.annotation.gtf.gz'),
#             'genomeInfo_file' = test_path('../../test_folders/static-data/UCSC_hg38_chromosome-info.tsv')
#         ),
#         'settings' = list(
#             'CNV_processing' = list(
#                 'gene_overlap' = list(
#                     'exclude_gene_type_regex' = c(),
#                     'include_only_these_gene_types' = c('lncRNA', 'miRNA', 'protein_coding'),
#                     'high_impact_list' = '__inbuilt__/supplemental-files/HighImpact-stemcell-hotspots.tsv',
#                     'highlight_list' = '__inbuilt__/supplemental-files/genelist-cancer-drivers.tsv'
#                 )
#             )
#         )
#     )
# 
#     high_impact_tb <- load_hotspot_table(config)
#     highlight_tb <- load_hotspot_table(config, 'Highlight')
#     config$settings$CNV_processing$gene_overlap$highlight_list <- '__inbuilt__/supplemental-files/genelist-cancer-hotspots.tsv'
#     highlight_tb2 <- load_hotspot_table(config, 'Highlight')
#     gr_info <- load_genomeInfo(config)
#     gr_genes <- load_gtf_data(config)
# 
#     parse_hotspot_table(high_impact_tb, gr_genes, gr_info)
#     parse_hotspot_table(highlight_tb, gr_genes, gr_info)
#     parse_hotspot_table(highlight_tb2, gr_genes, gr_info)
# 
# })


# 1 - not hit
# 2 - gene hit, any (but not gband): DDX11L1
# 3 - gene hit, gain: dummyC
# 4 - gene hit, no match
# 5 - position hit (any) + gband_hit (loss): chr1:40000-50000 + 1p36
# 6 - gband hit (any subband): 1p36
# 7 - gband hit (specific subband): 1p35.2

sample_cnvs <- tibble(
    seqnames = 'chr1',
    start = c(4000, 12000, 28050000, 28060000, 45000, 5000000, 31000000),
    end   = c(5500, 14000, 28055000, 28065000, 50000, 7000000, 32000000),
    sample_id = 'test_sample',
    CNV_type = c('gain', 'gain', 'gain', 'loss', 'loss', 'loss', 'loss'),
    ID = paste('combined', CNV_type, seqnames, start, end, sep='_'),
    CNV_caller = c('StemCNV-check', 'toolA', 'StemCNV-check', 'StemCNV-check', 'toolA', 'toolB', 'toolB'),
    n_probes = c(15, 5, 15, 25,5, 5, 5),
    CN = c(3, 3, 4, 1, 1, 1, 0),
    reference_overlap = c(T, T, F, F, F, F, T),
    reference_coverage = c(100, 85.01, NA_real_, NA_real_, NA_real_, NA_real_, 60),
    reference_caller = c('StemCNV-check', 'faketool', NA_character_, NA_character_,NA_character_,NA_character_, 'toolA')
) %>% as_granges()


#annotate_impact_lists
test_that('annotate_impact_lists', {
    gr_info <- load_genomeInfo(config)
    gr_genes <- load_gtf_data(config)
    
    hotspots <- parse_hotspot_table(read_tsv(test_path('../data/minimal-hotspots.tsv')), gr_genes, gr_info)
    
    expected_gr <- sample_cnvs %>%
        mutate(test_hits = c(NA, 'DDX11L1', 'dummyC', NA, '1p36|chr1:40000-50000', '1p36', '1p35.2'))
    expect_equal(annotate_impact_lists(sample_cnvs, hotspots, 'test'), expected_gr)

    # test empty hotspots
    expected_gr <- sample_cnvs %>%
        mutate(test_hits = NA_character_)
    expect_equal(annotate_impact_lists(sample_cnvs, GRanges(), 'test'), expected_gr)

})


test_that('annotate_roi', {
    gr_info <- load_genomeInfo(config)
    gr_genes <- load_gtf_data(config)
    
    sample_tb <- tibble(
        Sample_ID = 'test_sample',
        Regions_of_Interest = 'chr1:10000-20000;nametest|chr1:5000000-5500000;chrX:30000-40000;1p35'
    )

    expected_gr <- sample_cnvs %>%
        mutate(ROI_hits = c(NA, 'ROI_1', 'ROI_4', 'ROI_4', NA, 'nametest', 'ROI_4'))

    annotate_roi(sample_cnvs, 'test_sample', sample_tb, gr_genes, gr_info)
    
    expect_equal(annotate_roi(sample_cnvs, 'test_sample', sample_tb, gr_genes, gr_info), expected_gr)
    
    # test empty roi
    expected_gr <- sample_cnvs %>%
        mutate(ROI_hits = NA_character_)
    sample_tb$Regions_of_Interest <- ''
    expect_equal(annotate_roi(sample_cnvs, 'test_sample', sample_tb, gr_genes, gr_info), expected_gr)
    sample_tb <- tibble(Sample_ID = 'test_sample')
    expect_equal(annotate_roi(sample_cnvs, 'test_sample', sample_tb, gr_genes, gr_info), expected_gr)
})
