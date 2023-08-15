#! /usr/bin/Rscript
# Filter SNP probes based on settings in config
suppressMessages(library(argparse))

parser <- ArgumentParser(description="Filter SNP probes by quality scores")

parser$add_argument('inputfile', type = 'character', help='Path to input file')
parser$add_argument('outputfile', type = 'character', help='Path to output file')
parser$add_argument('configfile', type = 'character', help='Path to config file')

parser$add_argument('-f', '--filter-set', type = 'character', default = 'extended',
					help="Which filter-set (defined in config.yaml) to use")


args <- parser$parse_args(c('test/data/BIHi005-A/BIHi005-A.processed-data.tsv', 'test/filter-test.out.tsv', 'default_config.yaml'))
args <- parser$parse_args()

suppressMessages(library(tidyverse))
suppressMessages(library(yaml))

inputfile <- args$inputfile
outputfile <- args$outputfile
config <- read_yaml(args$configfile)
filter.set <- args$filter_set

if (!filter.set %in% names(config$settings$`probe-filter-sets`)) {
	stop(str_glue('Error: the probe-filter-set {filter.set} is not defined in the given config file.'))
}

GT.th        <- config$settings$`probe-filter-sets`[[filter.set]][['GenTrainScore']]
if (!is.numeric(GT.th) | GT.th >= 1 | GT.th <= 0) {
	warning(str_glue('threshold for GenTrainScore in filter-set "{filter.set}" is not numeric or out of bounds (0,1). Ignoring this filter.'))
	GT.th <- 0
}

GC.th        <- config$settings$`probe-filter-sets`[[filter.set]][['GenCallScore']]
if (!is.numeric(GC.th) | GC.th >= 1 | GC.th <= 0) {
	warning(str_glue('threshold for GenCallScore in filter-set "{filter.set}" is not numeric or out of bounds (0,1). Ignoring this filter.'))
	GC.th <- 0
}

multi.probes <- config$settings$`probe-filter-sets`[[filter.set]][['Position.duplicates']]
if (!multi.probes %in% c('keep', 'remove', 'highestGT', 'highestGC')) {
	warning(str_glue('value for Position.duplicates in filter-set "{filter.set}" is not onw of the allowed settings (keep|remove|highestGC|highestGT). Ignoring this filter.'))
	multi.probes.th <- 'keep'
}


tb <- read_tsv(inputfile, show_col_types = FALSE) %>%
	filter(`GenTrain Score` > GT.th & if_any(ends_with('gencall.Score'), ~ . > GC.th))

if (multi.probes != 'keep') {
	tb <- tb %>% group_by(Chr, Position) %>%
		mutate(multi = n() > 1)

	if (multi.probes == 'remove') {
		tb <- tb %>%
			ungroup() %>%
			filter(!multi) %>%
			dplyr::select(-multi)
	} else {
		select_by_col <- ifelse(multi.probes == 'highestGT', 'GenTrain Score',
		                        colnames(tb) %>% str_subset('gencall.Score$'))
		tb <- tb %>%
			slice_max(!!sym(select_by_col)) %>%
			dplyr::select(-multi)
	}
}

write_tsv(tb, outputfile)