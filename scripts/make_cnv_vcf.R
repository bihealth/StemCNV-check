#! /usr/bin/Rscript
# Convert tabular cnv info into vcf format
suppressMessages(library(tidyverse))
suppressMessages(library(plyranges))
suppressMessages(library(argparse))
suppressMessages(library(vcfR))
suppressMessages(library(yaml))

parser <- ArgumentParser(description="Convert tsv file with CNV calls from process_CNV_calls to vcf")

parser$add_argument('data_path', type = 'character', help='Path in which folders with sample data are located')
parser$add_argument('sample_id', type = 'character', help='CNV-pipeline sample_id')
parser$add_argument('config_path', type = 'character', help='Path to config file')

parser$add_argument('-m', '--mode', type = 'character', default = 'combined-calls',
					choices = c('split-tools', 'combined-calls'),
					help="Create VCF file for combined calls (default) or one per tool used.")
parser$add_argument('-i', '--include-state', type = 'character', nargs = '+',
					default = c('loss', 'gain', 'LOH'), choices = c('loss', 'gain', 'LOH'),
					help="Types of CNVs to include in VCF: gain|loss|LOH (default: all)")

args <- parser$parse_args()
#args <- parser$parse_args(c('data', '206210670080_R09C02', 'default_config.yaml', '-m', 'split-tools', '-i', 'gain', 'loss'))


sample_id <- args$sample_id
config <- read_yaml(args$config_path)
use.filter <- config$settings$filter$`use-filterset`
data_path <- args$data_path

# necessary info
# vcf contig headers -> can be taken from SNP vcf file (that way the genome model should be correct)
snp.vcf <- file.path(args$data_path, sample_id, paste0(sample_id, '.unprocessed.vcf')) %>%
	read.vcfR(verbose = F)

vcf.header <- snp.vcf@meta %>% str_subset('^##(fileformat|FILTER|contig)')

probe.data <- makeGRangesFromDataFrame(as_tibble(snp.vcf@fix) %>% select(-INFO, -ALT, -QUAL, -FILTER),
									   seqnames.field = 'CHROM', start.field = 'POS', end.field = 'POS',
									   keep.extra.columns = T, ignore.strand = T)

get_REF_entry <- function(seqname, pos) {
	gr <- filter_by_overlaps(probe.data,
							  GRanges(seqname, IRanges(start = pos, width = 1))) %>%
		as_tibble()
	if (nrow(gr) == 1) {
		return(gr$REF)
	} else if(nrow(gr) > 1) {
		warning(str_glue('Multiple Probes at CNV start: {seqname}:{pos}. Estimating REF allele from majority (alphabetic for ties).'))
		return(names(which.max(table(gr$REF))))
	} else {
		warning(str_glue('No Probes at CNV start: {seqname}:{pos}. This should not happen, can not determine REF allele.'))
		return('.')
	}
}



processed.calls <- file.path(data_path, sample_id,
							 paste0(sample_id, '.combined-cnv-calls.', use.filter, '.tsv')) %>%
	read_tsv(show_col_types = FALSE) %>%
	# Maybe needed for the list col redoing
	rowwise() %>%
	#Function to Get REF allele is not fast, but works
	mutate(
		#Need to reduce listcols somehow
		#TODO:
		# Max is not actually correct, but using sum here will be less accurate, this is not solvalble with the currently avalable information
		# Accurate number would need to be derived in parallel merging, requires loading the SNP subsets & doing a set merge on them
		numsnp = numsnp %>% str_split(';') %>% unlist() %>% max(),
		tool = tool %>% str_split(';') %>% unlist() %>%unique() %>% sort %>% paste(collapse='&'),
		# Only difference here will by 3/4/... or 0/1
		# -> should maybe take the PennCNV one (instead of majority / 3 or 0)?
		# -> Actually 3 or 0 are better assumptions if both 3/4 and 0/1 are the options.
		copynumber = copynumbers %>% str_split(';') %>% unlist() %>%table() %>% which.max %>% names(),
		CNV.type = ifelse(CNV.state == 'gain', 'DUP', 'DEL'),
		CNV.type = ifelse(CNV.state == 'LOH', 'LOH', CNV.type),
		CHROM = factor(seqnames, levels = paste0('chr', c(1:22, 'X', 'Y'))),
		POS = start, #TODO not sure if this is also 1-indexed
		ID = str_glue("CNV_{seqnames}_{start}_{end}"),
		REF = map2_chr(seqnames, start, get_REF_entry),
		ALT = str_glue('<{CNV.type}>'),
	    QUAL = '.',
		FILTER = '.',
		INFO = str_glue("END={end};SVTYPE={CNV.type};SVLEN={length}"),
		genotype = ifelse(CNV.state == 'LOH', '0/1', '1/1')
	 )

# MAYBE:
# If the processed table would contain start (&end) Probe of a CNV, a faster lookup than `filter_by_overlaps` could be used
# -> PennCNV would give that info as is; CBS not, can't know for future tools
# -> Therefore would require doing this probe matching step in the CBS &/ CNV processing steps
# -> Reading SNP probe data in CNV processing might allow easier merging by probe distance
# Overall:
# > need to load vcf here anyway since the start Base is irrelevant for the cnv table
# > might move the slow lookup to another script when including probe based merging into cnv-processing


#Static data for vcf writing
info.lines <- c('##INFO=<ID=END,Number=1,Type=Integer,Description="End coordinate of the variant">',
				'##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">',
				'##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">')

use.cols <- list(c('genotype', 'GT', "Segment genotype"),
				 c('copynumber', 'CN', "Segment most-likely or estimated copy-number call"),
				 c('numsnp', 'PN', "Number of points (i.e. SNP probes) in the segment"),
				 c('tool', 'TO', "Segment called as CNV by these tools")
)
vcf_cols <- c('CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT')

write_to_vcf <- function(tb, outvcf) {

	format.lines <- sapply(use.cols, function (x) {
		type <- ifelse(is.numeric(tb[[x[1]]]), 'Integer', 'String')
		shorthand <- x[2]
		desc <- x[3]
		str_glue('##FORMAT=<ID={shorthand},NUMBER=1,Type={type},Description="{desc}">')
	})

	writeLines(c(vcf.header, info.lines, format.lines,
	             paste0('#', paste(c(vcf_cols, sample_id), collapse = '\t'))), con = outvcf)


	tb %>%
		filter(CNV.state %in% args$include_state) %>%
		rowwise() %>%
		mutate(
			FORMAT = paste(sapply(use.cols, function (x) x[2]), collapse = ':'),
			#TODO this should be possibl dynamically somehow
			#sample_formatted = paste(sapply(use.cols, function (x) x[1]), collapse = ':'),
			sample_formatted = paste(genotype, copynumber, numsnp, tool, sep = ':'),
		) %>%
		dplyr::select(one_of(c(vcf_cols, 'sample_formatted', 'sample_id'))) %>%
		pivot_wider(names_from = sample_id, values_from = sample_formatted) %>%
		arrange(CHROM, POS) %>%
		write_tsv(outvcf, append=T)

	}


if (args$mode == 'split-tools') {
	lapply(config$settings$CNV.calling.tools, function(use.tool) {
		outvcf <- str_glue('{data_path}/{sample_id}/{sample_id}.{use.tool}-cnv-calls.{use.filter}.vcf')
		processed.calls %>%
			filter(tool.overlap.state != 'post-overlap' & tool == use.tool) %>%
			write_to_vcf(., outvcf = outvcf)
	})
} else {
	#TODO This should work, but there are some issues with duplicted calls? might only be GADA, unclear
	outvcf <- str_glue('{data_path}/{sample_id}/{sample_id}.combined-cnv-calls.{use.filter}.vcf')
	processed.calls %>%
		filter(tool.overlap.state != 'pre-overlap') %>%
		#some purrr walk function instead?
		write_to_vcf(., outvcf = outvcf)
}