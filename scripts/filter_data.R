#! /usr/bin/Rscript
# Filter SNP probes by quality scores
suppressMessages(library(argparse))

parser <- ArgumentParser(description="Filter SNP probes by quality scores")

parser$add_argument('inputfile', type = 'character', help='Path to input file')
parser$add_argument('outputfile', type = 'character', help='Path to output file')

parser$add_argument('-f', '--filter-setting', type = 'character', default = 'full',
					choices = c('basic', 'highGT', 'highGC', 'full'),
					help="Value for split SD undo")

args <- parser$parse_args()

suppressMessages(library(tidyverse))

inputfile <- args$inputfile
outputfile <- args$outputfile
filter.settings <- args$filter_setting


filters <- list(
	basic = function(x) dplyr::filter(x, across(6, ~ . > 0.15)),
	highGT = function(x) dplyr::filter(x, across(6, ~ . > 0.15) & across(12, ~ . > 0.8)),
	highGC = function(x) dplyr::filter(x, across(6, ~ . > 0.8)),
	full = function(x) dplyr::filter(x, across(6, ~ . > 0.8) & across(12, ~ . > 0.8))
)

ffunc <- filters[[filter.settings]]

read_tsv(inputfile, show_col_types = FALSE) %>%
	ffunc() %>%
	#TODO add splitting of the XY probes ?
	# 	mutate(sample_id = str_remove(filename, ....),
	# 			 Chr = paste0('chr', Chr),
	# 			 Chr = ifelse(Chr == 'chrXY', 'chrX,chrY', Chr)) %>%
	# separate_rows(Chr, sep=',')
	write_tsv(outputfile)
