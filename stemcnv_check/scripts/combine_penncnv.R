# Redirect all output to snakemake logging
log <- file(snakemake@log[['err']], 'wt')
sink(log, append = T)
sink(log, append = T, type = 'message')

library(tidyverse)
library(plyranges)
library(vcfR)
library(GenomeInfoDb)

snakemake@source('R/vcf_io_functions.R')
snakemake@source('R/processCNVs_calls_prefilter_merge.R')


read_PennCNV <- function(filename, sample_id) {
    read.table(filename, sep='', header = F, fill=T,
                         col.names = c('Position', 'numsnp', 'length', 'hmm.state', 'input', 'startsnp', 'endsnp', 'caller_confidence')) %>%
        separate(Position, c('seqnames', 'start', 'end'), convert=T) %>%
        # Note: start/end from PennCNV is from vcf, so also 1-based; length is start&end inclusive like granges
        # Note: this will go through granges before goint into vcf
        # dplyr::rename(sample_id = input, n_snp_probes = numsnp) %>%
        mutate(
            across(c(4,5,8,9,10), ~ str_remove(., '.*=')),
            across(c(4,10), ~as.numeric(.)),
            width = str_remove_all(length, ',') %>% as.integer(),
            sample_id = sample_id, #str_remove(input, '.*/') %>% str_remove('\\.filtered-data-.*\\.tsv$'),
            CNV_caller = 'PennCNV',   
            #TODO: this might be wrong for Y-chrom calls!! (or male X)
            CN = str_extract(hmm.state, '(?<=cn=)[0-9]') %>% as.integer(),  
            CNV_type = ifelse(CN < 2, 'DEL', NA),
            CNV_type = ifelse(CN == 2, 'CNV:LOH', CNV_type),
            CNV_type = ifelse(CN > 2, 'DUP', CNV_type),
            CNV_type = as.character(CNV_type),                      
            ID = paste(CNV_caller, CNV_type, seqnames, start, end, sep='_'),
        ) %>%
        select(seqnames, start, end, width, ID, CNV_caller, CNV_type, CN, sample_id)
}

penncnv_calls_to_vcf <- function (input_files, config, sample_id = 'test') { 
    # Config settings
    tool_config <- config$settings$PennCNV
    merge_distance <- tool_config$merge.distance

    #TODO: most of this can be made re-usable for CBS
    #TODO: code could be shortened a bit by functionaloising things
    
    #get SNP input vcf (with filters)
    snp_vcf <- read.vcfR(input_files[['vcf']], verbose = F) 
    snp_vcf_meta <- snp_vcf@meta
    # VCF POS should be 1-based,
    # Granges are also 1-based
    # AND are (fully) open [= start & end are included]
    snp_vcf_gr <- snp_vcf %>% 
        # only load necessary vcf cols, even if we get little to no speed-up
        vcfR_to_tibble(info_fields = F, format_fields = c("LRR", "BAF")) %>%
        as_granges(seqnames = CHROM, start = POS, width = 1) %>%
        filter(FILTER == 'PASS')
    
    # Read & pre-process all input files
    # Note: no more option to filter first & merge
    all_calls <- lapply(input_files[['tsvs']], read_PennCNV, sample_id = sample_id) %>%
        bind_rows() %>%
        as_granges()
    
    #Note: While PennCNV seems to require NCBI style CHROM in PFB file and input
    # it will still write out UCSC style in the ouput
    # -> need to convert from UCSC style to whatever style is used in SNP vcf
    styles <- genomeStyles() %>% .$Homo_sapiens
    snp_vcf_style <- seqlevelsStyle(snp_vcf_gr) %>% head(1)
    if (config$settings$vcf_output$chrom_style == '__keep__') {
        target_style <- snp_vcf_style     
    } else {
        target_style <- config$settings$vcf_output$chrom_style
    }
    # If target style is not alredy used in SNP vcf, need to fix it
    if (snp_vcf_style != target_style) {
        seqlevels(snp_vcf_gr) <- mapSeqlevels(
            seqnames(snp_vcf_gr) %>% .@values %>% as.character(), 
            target_style
        ) %>% sortSeqlevels()
    }   
    # Ensure that the fixed CHROM style from PennCNV is converted to target style
    if (seqlevelsStyle(all_calls) %>% head(1) != target_style) {
        seqlevels(all_calls) <- mapSeqlevels(
            seqnames(all_calls) %>% .@values %>% as.character(), 
            target_style
        ) %>% sortSeqlevels() %>% unique()
    }
        
    all_calls <- all_calls %>%
        # Properly sort
        arrange(seqnames, start) %>%
        # Combined nearby calls & update SNP counts
        merge_calls(merge_distance, snp_vcf_gr) %>%
        # Add FILTER column
        add_call_prefilters(tool_config) %>%
        arrange(seqnames, start) %>%
        as_tibble()        
    
    filtersettings <- tool_config$`filter-settings`
    if (filtersettings == '__default__') {
        filtersettings <- config$settings$`default-filter-settings`
    }    
    enable_LOH_calls <- tool_config$enable_LOH_calls
    header <- c(
        str_subset(snp_vcf_meta, 'fileformat|contig|BPM=|EGT=|CSV='),
        static_cnv_vcf_header(tool_config),
        #Add a line describing tool specific details
        '##ALT=<ID=CNV:LOH,Description="Loss of heterozygosity, same as run of homozygosity">',
        paste(
            '##PennCNV="docker://genomicslab/penncnv"',
            str_glue('StemCNV-check_array_probe_filtering="{filtersettings}"'),
            ifelse(enable_LOH_calls, 'LOH_detection=True', 'LOH_detection=False'),
            'GC_wave_correction=True'
        )
    )
    
    cnv_vcf <- new(
        "vcfR",
        meta = header,
        fix = get_fix_section(all_calls),
        gt = get_gt_section(all_calls, snp_vcf_gr)
    )
    
    write.vcf(cnv_vcf, snakemake@output[['vcf']])
  
}

penncnv_calls_to_vcf(snakemake@input, snakemake@config, snakemake@wildcards[['sample_id']])