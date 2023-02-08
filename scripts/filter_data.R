suppressMessages(library(tidyverse))
suppressMessages(library(optparse))

parser <- OptionParser(
	usage = "usage: %prog /path/to/inputfile.tsv /path/to/outputfile.tsv"
)

parser <- add_option(parser, c("-f", "--filter-setting"), 
					 default='full', help="Chose filter settings")

args <- parse_args(parser, positional_arguments = 2)

inputfile <- args$args[1]
outputfile <- args$args[2]
filter.settings <- args$options$`filter-setting`


filters <- list(
	basic = function(x) dplyr::filter(x, across(6, ~ . > 0.15)),
	highGT = function(x) dplyr::filter(x, across(6, ~ . > 0.15) & across(12, ~ . > 0.8)),
	highGC = function(x) dplyr::filter(x, across(6, ~ . > 0.8)),
	full = function(x) dplyr::filter(x, across(6, ~ . > 0.8) & across(12, ~ . > 0.8))
)

ffunc <- filters[[filter.settings]]

read_tsv(inputfile, show_col_types = FALSE) %>%
	ffunc() %>%
	write_tsv(outputfile)
