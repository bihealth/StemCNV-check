suppressMessages(library(tidyverse))
suppressMessages(library(plyranges))
suppressMessages(library(yaml))

defined_labels <- yaml.load_file(test_path('../../src/stemcnv_check/control_files/label_name_definitions.yaml'))

# Example calls (toolA and toolB):

# 0+1 - unqiue in one tool
# 1+0 - unqiue in one tool
# 1+1 - same call in both tools
# 1+1 - same region, but diff CNV_type
# 1+1 - both tools, but median region cov <60% (100%: 2k, <20% : <400)
# 1+1 - same call, but different copynumber

# 2+1 - 3 calls, all merged
# 2+1 - 3 calls, region with 2 calls; median region cov <60%  (100%: 9k, <20% : <1.8k)
# 2+1 - 3 calls, region with 1 call ; single call <50% region cov (4k---4k + 3k)
# 2+2 - 4 calls, both regions >60% median cov, but no single call >50% cov

toolA <- tibble(
    seqnames = rep('chr1', 13),
    start = c(1000, 4000, 6000,  8000, 120000, 15e5, 18e5, 22.2e4, 27e4, 31e4, 36e4, 50e4, 55e4) %>% as.integer(),
    end   = c(2000, 5500, 7500, 10000, 140000, 17e5, 20e5, 23e4, 27.8e4, 35e4, 40e4, 54e4, 58e4) %>% as.integer(),
    CNV_type = c(rep('DUP', 3), rep('DEL', 10)),
    sample_id = 'test_sample',
    FILTER = c('min_size', 'Test;Test2', 'min_probes', rep(NA_character_, 10)),
    CNV_caller = 'toolA',
    n_probes = 10,
    n_uniq_probes = n_probes,
    probe_density_Mb = n_uniq_probes / (end - start + 1) * 1e6,
    n_initial_calls = 1,
    initial_call_details = NA_character_,
    CN = c(rep(3, 3), rep(1, 10)),
    ID = paste(CNV_caller, CNV_type, seqnames, start, end, sep='_'),
)

toolB <- tibble(
    seqnames = rep('chr1', 10),
    start = c(2000, 4000, 6000,  9700, 120000, 15e5, 21e4, 33e4, 52e4, 57e4) %>% as.integer(),
    end   = c(3000, 5500, 7500, 10000, 140000, 20e5, 30e4, 36e4, 56e4, 60e4) %>% as.integer(),
    CNV_type = c(rep('DUP', 2), rep('DEL', 8)),
    sample_id = 'test_sample',
    FILTER = c(NA_character_, NA_character_, 'min_density', rep(NA_character_, 7)),
    CNV_caller = 'toolB',
    n_probes = 5,
    n_uniq_probes = n_probes,
    probe_density_Mb = n_uniq_probes / (end - start + 1) * 1e6,
    n_initial_calls = 1,
    initial_call_details = NA_character_,
    CN = c(rep(3, 2), 1, 1, 0, rep(1, 5)),
    ID = paste(CNV_caller, CNV_type, seqnames, start, end, sep='_'),
)


combined_tools <- bind_ranges(
    #combined calls
    # sortebd by: CNV_type (DEL before DUP), start
    tibble(     
        seqnames = 'chr1',
        start = c(120000, 15e5, 4000) %>% as.integer(),
        end   = c(140000, 20e5, 5500) %>% as.integer(),
        sample_id = 'test_sample',
        CNV_type = c('DEL', 'DEL', 'DUP'),
        ID = paste(defined_labels$combined_cnvs, CNV_type, seqnames, start, end, sep='_'),
        n_initial_calls = c(2, 3, 2),
        initial_call_details = c(
            "toolA_120000-140000_CN1_cov100_PASS|toolB_120000-140000_CN0_cov100_PASS",
            "toolA_1500000-1700000_CN1_cov40_PASS|toolA_1800000-2000000_CN1_cov40_PASS|toolB_1500000-2000000_CN1_cov100_PASS",
            "toolA_4000-5500_CN3_cov100_Test&Test2|toolB_4000-5500_CN3_cov100_PASS"
        ),       
        CNV_caller = defined_labels$combined_cnvs,
        CN = c(0.5, 1, 3),
        # overlap_merged_call = NA_real_,
        # Recalculated based on vcf file; indiv call will keep whatever (fake) number they had before
        n_probes = c(11, 6, 5),
        n_uniq_probes = n_probes,
        # this *should* be correct
        probe_density_Mb = n_uniq_probes / (end - start + 1) * 1e6,
        FILTER = c(NA, NA, 'min_probes'),
    ) %>% as_granges(),
    # single calls
    # sorted by type, then caller, then start
    bind_rows(
        toolA[c(1,3,4,8:13),],
        toolB[c(1,3,4,7:10),]
    ) %>%
        arrange(CNV_type, CNV_caller, start) %>%
        mutate(
                n_initial_calls = 1,
                initial_call_details = NA_character_,
        ) %>%
        as_granges()
)