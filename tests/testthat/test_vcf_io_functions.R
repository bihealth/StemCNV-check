library(testthat)

library(tidyverse)
library(plyranges)
library(GenomeInfoDb)
library(vcfR)

source(test_path('../../src/stemcnv_check/scripts/R/helper_functions.R'))
source(test_path('../../src/stemcnv_check/scripts/R/vcf_io_functions.R'))

#  Functions to test:
# - [x] vcfR_to_tibble
# - [x] parse_snp_vcf
# - [x] parse_cnv_vcf
# - [x] get_fix_section
# - [x] get_gt_section
# - [x] write_cnv_vcf (incl. vcf header)
# - [x] fix_header_lines

snp_vcfr <- read.vcfR(test_path('../data/minimal_probes.vcf'), verbose = F) 

expected_snp_tb <- tibble(
    CHROM = 'chr1',
    POS = c(1, 49, 100, 105, 110, 115, 115, 199, 1000, 1300, 1599, 3000, 3250, 3500, 3500, 3750, 3999) %>% as.integer(),
    ID = c('dummy1a', 'dummy1b', NA, NA, NA, NA, 'dummy', NA, NA, NA, NA, NA, NA, NA, 'DUP', NA, NA),
    REF = c('A', 'C', 'T', 'G', 'A', 'C', 'C', 'G', 'A', 'C', 'T', 'G', 'A', 'T', 'T', 'G', 'C'),
    ALT = c('C', 'T', 'G', 'A', 'C', 'T', 'T', 'A', 'C', 'T', 'G', 'A', 'C', 'G', 'G', 'A', 'T'),
    QUAL = NA_real_, # or "."
    FILTER = c(rep('PASS', 6), 'LOWQUAL', rep('PASS', 10)),
    GenTrain_Score = c(rep(0.8, 6), 0.2, rep(0.8, 7), 0.9, 0.8, NA),
    CSQ = c(NA,NA,'gene_id|gene_name', NA, NA, NA, NA, NA, NA, NA, NA, 'ABC|abc', NA, NA, NA, NA, NA),
    sample_id = 'test_sample',
    GT = c(rep('0/1', 16), '0/0'),
    LRR = c(1,1,1, 1.1, 0.7, 1, 1, 0.9, 1.3, 1.2, 1.3, 0.85, 0.82, 0.8, 1.15, 0.87, 0.88),
    BAF = c(1,1,1, 0, 0.3, 0.7, 1, 1, 1, 1, 0, 0.33, 1, 1, 0.67, 1, 0.67)        
)

test_that("vcfR_to_tibble", {
    vcfR_to_tibble(snp_vcfr) %>%
        filter(POS < 4000) %>% 
        expect_equal(expected_snp_tb)
})

test_that('parse_snp_vcf', {

    expected <- expected_snp_tb %>%
        dplyr::rename(seqnames = CHROM, start = POS) %>%
        mutate(width = 1) %>%
        as_granges()
    
    parse_snp_vcf(snp_vcfr, apply_filter = FALSE, info_fields = NULL, format_fields = NULL) %>%
        filter(start < 4000) %>% 
        expect_equal(expected)    
    
    expected <- expected %>%
        select(-GenTrain_Score, -CSQ, -GT) %>%
        filter(FILTER == 'PASS')
    
    parse_snp_vcf(test_path('../data/minimal_probes.vcf')) %>%
        filter(start < 4000) %>% 
        expect_equal(expected)

})

cnv_tb <- tibble(
    #vcf/granges block
    seqnames  = c('chr1', 'chr1', 'chr1', 'chr3', 'chr5', 'chr17', 'chr18', 'chrX', 'chrX') %>%
        factor(levels = genomeStyles('Homo_sapiens')$UCSC),
    start = c(100, 1000, 3000, 6000, 1.0e8, 7400,  9000, 12e4, 1e6) %>% as.integer(),
    end   = c(200, 1600, 5400, 7000, 1.1e8, 8400, 10000, 15e4, 3e6+400)-1 %>% as.integer(),
    width = c(100, 600, 2400, 1000, 1e7, 1000, 1000, 3e4, 2e6+400),
    #ID goes here
    # REF & QUAL removed; CNV_type = ALT
    CNV_type = c(rep('gain', 5), 'LOH', rep('loss', 2), 'gain'),  
    
    FILTER = c('min_size', 'min_size;min_probes', NA, NA, 'min_density', NA, NA, NA, 'min_density'),
    # INFO block
    n_probes =      c(5, 3, 11, 5, 20, 5, 10, 15, 20),
    n_uniq_probes = c(5, 3, 10, 5, 20, 5, 10, 15, 20),
    probe_density_Mb = round(n_uniq_probes / width * 1e6, 3),
    # Check_Score, Call_label, stemcell_hotspot, dosage_sensitive_gene, cancer_gene, Gap_percent, overlapping_genes
    # FORMAT & sample block  
    sample_id = 'test_sample',
    # GT is removed
    CN = c(3, 4, 3, 3, 3, 2, 0, 1, 3),  
    LRR = c(1, 1.3, 0.89, 2, 1.385, 0, -0.88, -1.31, -0.895),
    #(for the future:) BAF_clusters = c(4, 2, 4, 3, 4, 2, 2, 2, 2)
    CNV_caller = 'Test',  
    n_initial_calls = c(1, 1, 2, 1, 1, 1, 2, 3, 2),
    initial_call_details = c(
        NA, NA, '3000_3999_CN3|4400_5399_CN3', NA, NA, NA, '9000_9399_CN2|9600_9999_CN2', 
        '120000_129999_CN1,130000_139999_CN0|140000_149999_CN1', '1000000_1999999_CN1|2000400_3000399_CN1'
    ),
    # ROI_hits, precision_estimate, precision_estimate_description, reference_coverage
    ID = paste(CNV_caller, CNV_type, seqnames, start, end, sep='_'),
) %>%
    relocate(ID, .after = width)

#Add additional info for fully processed/annotated CNV calls
cnv_tb_annotated <- cnv_tb %>%
    mutate(
        # extra INFO block
        # FILTER = c(), # could be adapted to proper values, but doesn't matter for this test
        Check_Score = c(68.924, 38.838, 45.692, 37.664, 25.518, 26.937, 30.998, 32.161, 47.375),
        Call_label = c('Critical', NA, 'Reportable', NA, NA, NA, NA, rep ('Reference genotype', 2)),
        stemcell_hotspot = c(NA, NA, NA, NA, NA, 'gene1,gene2', NA, NA, NA),
        dosage_sensitive_gene = NA_character_,
        cancer_gene = c(NA, NA, 'gene3', NA, NA, NA, NA, NA, NA),
        Gap_percent = c(0, 0, 0, c(0.943, 0.047, 0.102, 0.212, 0.303, 0.652)), # random numbers, hard-coded to vcf
        # not actually used in the function
        # high_probe_density = c(FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, TRUE),
        overlapping_genes = c('abc,er123,xyz', NA, 'gene3', NA, NA, 'gene1,gene2', NA, NA, 'gene5'),
        # extra FORMAT block
        ROI_hits = c(NA, NA, NA, NA, NA, NA, NA, NA, 'ROI1'),
        precision_estimate = c(0.6, 0.4, 0.4, 0.8, 0.6, 0.1, 0.4, NA, NA),
        precision_estimate_description = 'dummy: somevalue; another: value',
        # reference_caller
        reference_coverage = c(rep(NA, 4), c(0.239, 0.208, 0.246, 0.777, 0.812)), # random numbers, hard-coded to vcf
    ) %>%
    relocate((ncol(cnv_tb)+1):(ncol(cnv_tb)+7), .after = probe_density_Mb)

cnv_tb_annotated_out <- cnv_tb_annotated %>%
    mutate(
        across(where(is.numeric), ~ round(., 3) %>% as.character()),
        across(where(is.character), ~ ifelse(is.na(.), '.', .) %>% str_replace_all(',', '|')),
        precision_estimate_description = 'dummy=somevalue;another=value'
    )

cnv_tb_empty <- cnv_tb %>% filter(seqnames == 'dummy')

test_that('get_fix_section', {
    expected_fix <- tibble(
        CHROM = cnv_tb$seqnames,
        POS = cnv_tb$start - 1,
        ID = cnv_tb$ID,
        REF = '.',
        ALT = paste0(
            '<', cnv_tb$CNV_type %>% str_replace('gain', 'DUP') %>% 
                str_replace('loss', 'DEL') %>% str_replace('LOH', 'CNV:LOH'),
            '>'
        ),
        QUAL = '.',
        FILTER = ifelse(is.na(cnv_tb$FILTER), 'PASS', cnv_tb$FILTER),
        # 'END={end};SVLEN={width};SVCLAIM=D;N_PROBES={n_probes};N_UNIQ_PROBES={n_uniq_probes};PROBE_DENS={probe_density_Mb}'
        INFO = paste0(
            str_glue('END={cnv_tb$end};SVLEN={cnv_tb$width};'),
            'SVCLAIM=D;',
            str_glue('N_PROBES={cnv_tb$n_probes};N_UNIQ_PROBES={cnv_tb$n_uniq_probes};'),
            str_glue('PROBE_DENS={round(cnv_tb$probe_density_Mb, 3)}')            
        ),
    ) %>% 
        arrange(CHROM, POS) %>%
        mutate(POS = as.character(POS))
    
    get_fix_section(cnv_tb) %>%
        expect_equal(as.matrix(expected_fix))
    
    # test empty
    expected_empty <- expected_fix %>% 
        filter(CHROM == 'dummy') %>%
        as.matrix()
    
    get_fix_section(cnv_tb_empty) %>%
        expect_equal(expected_empty)
    
    # test advanced 
    expected_fix <- expected_fix %>%
        mutate(
            INFO = paste0(
                str_glue('END={cnv_tb$end};SVLEN={cnv_tb$width};'),
                'SVCLAIM=D;',
                str_glue('N_PROBES={cnv_tb$n_probes};N_UNIQ_PROBES={cnv_tb$n_uniq_probes};'),
                str_glue('PROBE_DENS={round(cnv_tb$probe_density_Mb, 3)};'),
                str_glue('Check_Score={cnv_tb_annotated_out$Check_Score};'),
                str_glue('Call_label={cnv_tb_annotated_out$Call_label};'),
                str_glue('stemcell_hotspot={cnv_tb_annotated_out$stemcell_hotspot};'),
                str_glue('dosage_sensitive_gene={cnv_tb_annotated_out$dosage_sensitive_gene};'),
                str_glue('cancer_gene={cnv_tb_annotated_out$cancer_gene};'),
                str_glue('Gap_percent={cnv_tb_annotated_out$Gap_percent};'),
                paste0('overlapping_genes=', cnv_tb_annotated_out$overlapping_genes %>% str_replace(',', '|'))
            )
        )
    
    expect_equal(get_fix_section(cnv_tb_annotated), as.matrix(expected_fix))
})


test_that('get_gt_section', {
    # FORMAT keys: GT, CN, LRR, TOOL (str desc)
    expected_gt <- tibble(
        FORMAT = 'GT:CN:LRR:TOOL',
        test_sample = paste(
            c('0/1', './.', '0/1', '0/1', '0/1', './.', '1/1', '0/1', '0/1'),
            cnv_tb$CN,
            cnv_tb$LRR,
            paste0(
                "caller=Test;",
                "n_initial_calls=", cnv_tb$n_initial_calls, ";",
                "initial_call_details=", ifelse(is.na(cnv_tb$initial_call_details), '.', cnv_tb$initial_call_details)
            ),            
            sep = ":"
        )
    ) 
    expected_empty <- expected_gt %>% 
        filter(FORMAT == 'dummy') %>%
        as.matrix()
    expected_gt <- expected_gt %>% as.matrix()
    
    get_gt_section(cnv_tb, 'test_sample', 'f', 'UCSC') %>%
        expect_equal(expected_gt)
    
    # test empty
    get_gt_section(cnv_tb_empty, 'test_sample', 'f', 'UCSC') %>%
        expect_equal(expected_empty)

    # Test with changes for male X & Y
    expected_gt_m <- expected_gt
    expected_gt_m[17] <- "./.:3:-1.31:caller=Test;n_initial_calls=3;initial_call_details=120000_129999_CN1,130000_139999_CN0|140000_149999_CN1"
    expected_gt_m[18] <- "1/1:0:-0.895:caller=Test;n_initial_calls=2;initial_call_details=1000000_1999999_CN1|2000400_3000399_CN1"
    cnv_tb %>%
        mutate(
            seqnames = c('chr1', 'chr1', 'chr1', 'chr3', 'chr5', 'chr17', 'chr18', 'chrX', 'chrY') %>%
                factor(levels = genomeStyles('Homo_sapiens')$UCSC),
            CN = c(3, 4, 3, 3, 3, 2, 0, 3, 0)
        ) %>%
        get_gt_section('test_sample', 'm', 'UCSC') %>%
        expect_equal(expected_gt_m)
    
    # test with processed annotation
    expected_gt[1:9] <- 'GT:CN:LRR:TOOL:ROI:PREC_EST:PREC_DESC:REFCOV'
    expected_gt[10:18] <- paste(
        expected_gt[10:18],
        cnv_tb_annotated_out$ROI_hits,
        cnv_tb_annotated_out$precision_estimate,
        cnv_tb_annotated_out$precision_estimate_description,
        cnv_tb_annotated_out$reference_coverage,
        sep = ':'
    )
    get_gt_section(cnv_tb_annotated, 'test_sample', 'f', 'UCSC') %>%
        expect_equal(expected_gt )
})

# fix_header_lines <- function(header_lines, regex = NULL, contig_format = NULL)
test_that('fix_header_lines', {
    expected <- c(
        "##contig=<ID=chr1,length=249250621>",
        "##contig=<ID=chr3,length=198022430>",
        "##contig=<ID=chr5,length=180915260>",
        "##contig=<ID=chr17,length=81195210>",
        "##contig=<ID=chr18,length=78077248>",
        "##contig=<ID=chrX,length=155270560>",
        "##BPM=manifest_A1.bpm",
        "##CSV=manifest_A1.csv.gz",
        "##EGT=manifest_A1.egt"
    )
    header_lines <- snp_vcfr@meta
    expect_equal(fix_header_lines(header_lines, 'contig|manifest', contig_format = 'UCSC'), expected)
})

#parse_cnv_vcf <- function(vcf, info_fields = NULL, format_fields = NULL, apply_filter = FALSE)
test_that('parse_cnv_vcf', {
    expected <- cnv_tb_annotated %>%
        mutate(
            across(
                c(stemcell_hotspot, dosage_sensitive_gene, cancer_gene, overlapping_genes),
                ~str_replace_all(., ',', '|')
            ),
            seqnames = factor(seqnames, levels = c('chr1', 'chr3', 'chr5', 'chr17', 'chr18', 'chrX'))
        ) %>%
        as_granges() 
    
    parse_cnv_vcf(test_path('../data/minimal-cnvs.vcf')) %>%
        expect_equal(expected)
})


# write_cnv_vcf <- function(cnv_tb, out_vcf, sample_sex, tool_name, fullconfig, defined_labels, snp_vcf_meta, command_desc, target_style) 
test_that('write_cnv_vcf', {
    fullconfig <- list(
        settings = list(
            default_probe_filter_set = 'standard',
            CNV_processing = list(
                call_processing = list(
                  probe_filter_settings = '_default_',
                  filter.minprobes = 5,
                  filter.minlength = 1000,
                  filter.mindensity.Mb = 10, #snps per Mb
                  min.perc.gap_area = 0.33,
                  density.quantile.cutoff = 0.99,
                  tool.overlap.greatest.call.min.perc = 50,
                  tool.overlap.min.cov.sum.perc = 60,
                  min.reciprocal.coverage.with.ref = 50
                )
            )
        )
    )

    tmp_vcf <- tempfile(fileext = '.vcf.gz')
    
    write_cnv_vcf(
        cnv_tb_annotated,
        tmp_vcf,
        'female',
        'combined_calls',
        fullconfig,
        defined_labels,
        snp_vcfr@meta,
        'StemCNV-check test description',
        'UCSC'
    )
    R.utils::gunzip(tmp_vcf)
    expect(
        compare_file_text(test_path('../data/minimal-cnvs.vcf'), str_remove(tmp_vcf, '\\.gz')),
        'vcf files are not equal'
    )
})

