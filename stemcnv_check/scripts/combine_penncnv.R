# Redirect all output to snakemake logging
log <- file(snakemake@log[['err']], 'wt')
sink(log, append = T)
sink(log, append = T, type = 'message')

library(tidyverse)
library(plyranges)
library(vcfR)

snakemake@source('R/vcf_io_functions.R')
snakemake@source('R/processCNVs_calls_prefilter_merge.R')

read_PennCNV <- function(filename) {
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
            sample_id = snakemake@wildcards[['sample_id']], #str_remove(input, '.*/') %>% str_remove('\\.filtered-data-.*\\.tsv$'),
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

penncnv_calls_to_vcf <- function (input_files, config) { 
    # Config settings
    tool_config <- config$settings$PennCNV
    merge_distance <- tool_config$merge.distance

    #get SNP input vcf (with filters)
    snp_vcf <- read.vcfR(input_files[['vcf']], verbose = F)
    #TODO: directly apply FILTER col?
    
    # Read & pre-process all input files
    # Note: no more option to filter first & merge
    all_calls <- lapply(input_files[['tsvs']], read_PennCNV) %>%
        bind_rows() %>%
        # Combined nearby calls & update SNP counts
        merge_calls(merge_distance, snp_vcf) %>%
        # Add FILTER column
        add_call_prefilters(tool_config) %>%
        as_tibble()
    
    filtersettings <- tool_config$`filter-settings`
    if (filtersettings == '__default__') {
        filtersettings <- config$settings$`default-filter-settings`
    }    
    enable_LOH_calls <- tool_config$enable_LOH_calls
    header <- c(
        str_subset( snp_vcf@meta, 'fileformat|contig|BPM=|EGT=|CSV='),
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
        gt = get_gt_section(all_calls, snp_vcf)
    )
    
    write.vcf(cnv_vcf, snakemake@output[['vcf']])
  
}

penncnv_calls_to_vcf(snakemake@input, snakemake@config)