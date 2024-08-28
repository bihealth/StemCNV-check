# Redirect all output to snakemake logging
log <- file(snakemake@log[['err']], 'wt')
sink(log, append = T)
sink(log, append = T, type = 'message')

library(tidyverse)
library(plyranges)
library(vcfR)
library(GenomeInfoDb)

snakemake@source('R/helper_functions.R')
snakemake@source('R/vcf_io_functions.R')
snakemake@source('R/preprocess_CNV_functions.R')


read_PennCNV <- function(filename, sample_id, sample_sex) {
    tb <- read.table(filename, sep='', header = F, fill=T,
                         col.names = c('Position', 'numsnp', 'length', 'hmm.state', 'input', 'startsnp', 'endsnp', 'caller_confidence')) %>%
        separate(Position, c('seqnames', 'start', 'end'), convert=T) %>%
        # Note: start/end from PennCNV is from vcf, so also 1-based; length is start&end inclusive like granges
        # Note: this will go through granges before goint into vcf
        # dplyr::rename(sample_id = input, n_snp_probes = numsnp) %>%
        mutate(
            across(c(4,5,8,9,10), ~ str_remove(., '.*=')),
            across(c(4,10), ~as.numeric(.)),
            seqnames = as.character(seqnames),
            width = str_remove_all(length, ',') %>% as.integer(),
            sample_id = sample_id, 
            CNV_caller = 'PennCNV',
        ) 

    # Checking seqlevelStyle won't work on empty values
    if (nrow(tb) == 0) {
        tb <- tb %>%
            mutate(CN = integer(), ID = character(), CNV_type = character()) %>%
            select(seqnames, start, end, width, ID, CNV_caller, CNV_type, CN, sample_id)
        
        return(tb)
    }
    
    sex_chroms <- get_sex_chroms(tb)

    tb %>% 
        mutate(
            CN = str_extract(hmm.state, '(?<=cn=)[0-9]') %>% as.integer(),  
            CNV_type = case_when(
                # Male sex chroms have different CN <> CNV_type relations
                sample_sex == 'm' & seqnames %in% sex_chroms & CN < 1  ~ 'DEL',
                sample_sex == 'm' & seqnames %in% sex_chroms & CN > 1  ~ 'DUP',
                sample_sex == 'm' & seqnames %in% sex_chroms & CN == 1 ~ 'ERROR',
                # autosomes & female X
                CN < 2                                                 ~ 'DEL',
                CN == 2                                                ~ 'CNV:LOH',
                #FIXME (future): technically DUP should ONLY be used for CN=3 (or 2 on male XY)
                CN > 2                                                 ~ 'DUP',
            ),
            ID = paste(CNV_caller, CNV_type, seqnames, start, end, sep='_')
        ) %>%
        select(seqnames, start, end, width, ID, CNV_caller, CNV_type, CN, sample_id)
}

penncnv_calls_to_vcf <- function(input_files, out_vcf, config, sample_id = 'test') { 
    # Config settings
    tool_config <- config$settings$PennCNV
    
    #get SNP input vcf (with filters)
    snp_vcf <- read.vcfR(input_files[['vcf']], verbose = F) 
    snp_vcf_meta <- snp_vcf@meta
    snp_vcf_gr <- parse_snp_vcf(snp_vcf)
    
    # Read input files
    sample_sex <- get_sample_info(sample_id, "sex", config$sample_table)
    all_calls <- lapply(input_files[['tsvs']], read_PennCNV, 
                        sample_id = sample_id, sample_sex = sample_sex) %>%
        bind_rows() %>%
        as_granges()
    
    #make sure that snp_vcf & cnv_vcf use the same (& intended) chrom style
    target_style <- get_target_chrom_style(config, snp_vcf_gr)
    all_calls <- fix_CHROM_format(all_calls, target_style)
    snp_vcf_gr <- fix_CHROM_format(snp_vcf_gr, target_style)
    
    # preprocess (merge, filter, SNP counts)
    all_calls <- apply_preprocessing(all_calls, snp_vcf_gr, tool_config) %>%
        get_median_LRR(snp_vcf_gr) %>%
        as_tibble()       
    
    filtersettings <- tool_config$`filter-settings`
    if (filtersettings == '__default__') {
        filtersettings <- config$settings$`default-filter-settings`
    }    
    enable_LOH_calls <- tool_config$enable_LOH_calls
    header <- c(
        fix_header_lines(snp_vcf_meta, 'fileformat|contig|BPM=|EGT=|CSV=', target_style),
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
        gt = get_gt_section(all_calls, sample_sex)
    )
    
    write.vcf(cnv_vcf, out_vcf)
  
}

penncnv_calls_to_vcf(
    snakemake@input,
    snakemake@output[['vcf']],
    snakemake@config,
    snakemake@wildcards[['sample_id']]
)