#! /usr/bin/Rscript
# Convert tabular cnv info into vcf format
suppressMessages(library(argparse))

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

suppressMessages(library(tidyverse))
suppressMessages(library(plyranges))
suppressMessages(library(vcfR))
suppressMessages(library(yaml))

sample_id <- args$sample_id
config <- read_yaml(args$config_path)
use.filter <- ifelse(config$settings$make_cnv_vcf$`filter-settings` == '__default__',
					 config$settings$`default-filter-set`,
					 config$settings$make_cnv_vcf$`filter-settings`)
name_addition <- config$settings$make_cnv_vcf$name_addition
data_path <- args$data_path

# necessary info
# vcf contig headers -> can be taken from SNP vcf file (that way the genome model should be correct)
snp.vcf <- file.path(args$data_path, sample_id, paste0(sample_id, '.unprocessed.vcf')) %>%
	read.vcfR(verbose = F)

vcf.header <- snp.vcf@meta %>% str_subset('^##(fileformat|FILTER|contig)')

probe.data <- makeGRangesFromDataFrame(as_tibble(snp.vcf@fix) %>% select(-INFO, -ALT, -QUAL, -FILTER),
									   seqnames.field = 'CHROM', start.field = 'POS', end.field = 'POS',
									   keep.extra.columns = T, ignore.strand = T)

# CHROM comes from GRanges and is alwayas 'chr1..22' -> need to adapt to format used in vcf header
remove_str <- ifelse(any(str_detect(vcf.header, '##contig=<ID=[0-9],')), 'chr', '---')
chrom_levels <- paste0(ifelse(remove_str == 'chr', '', 'chr'), c(1:22, 'X', 'Y'))

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
							 paste0(sample_id, '.combined-cnv-calls.tsv')) %>%
	read_tsv(show_col_types = FALSE) %>%
	unique() %>%
	# Maybe needed for the list col redoing
	rowwise() %>%
	#Function to Get REF allele is not fast, but works
	mutate(
		seqnames = str_remove(seqnames, remove_str),
		## Need to reduce listcols to a single value for vcf
		# Use highest tool_confidence  of single call
		caller_confidence = suppressWarnings(caller_confidence %>% str_split(';') %>% unlist() %>% as.numeric() %>% max()),
		#Note:
		# Max is not actually correct, but using sum here will be less accurate, this is not solvalble with the currently avalable information
		# Accurate number would need to be derived parallel to merging, requires loading the SNP subsets & doing a set merge on them
		n_snp_probes = n_snp_probes %>% str_split(';') %>% unlist() %>% as.integer() %>% max(),
		CNV_caller = CNV_caller %>% str_split(';') %>% unlist() %>% unique() %>% sort %>% paste(collapse='&'),
		# No good way to decied if multiple copynumbers are predicted
		# Only difference here will by 3/4/... or 0/1
		# -> should maybe take the PennCNV one (instead of majority / 3 or 0)?
		# -> Maybe 3 or 0 are better assumptions if both 3/4 and 0/1 are the options.
		copynumber = copynumber %>% str_split(';') %>% unlist() %>%table() %>% which.max %>% names(),
		CNV_type = ifelse(CNV_type == 'gain', 'DUP', 'DEL'),
		CNV_type = ifelse(CNV_type == 'LOH', 'LOH', CNV_type),
		CHROM = factor(seqnames, levels = chrom_levels),
		POS = start - 1, # 1-indexed, but:
		#https://samtools.github.io/hts-specs/VCFv4.2.pdf:
		# For simple insertions and
		# deletions in which either the REF or one of the ALT alleles would otherwise be null/empty, the REF and ALT
		# Strings must include the base before the event (which must be reflected in the POS field), unless the event
		# occurs at position 1 on the contig in which case it must include the base after the event; this padding base is
		# not required (although it is permitted) for e.g. complex substitutions or other events where all alleles have at
		# least one base represented in their Strings. If any of the ALT alleles is a symbolic allele (an angle-bracketed
		# ID String “<ID>”) then the padding base is required and POS denotes the coordinate of the base preceding
		# the polymorphism.
		ID = str_glue("CNV_{seqnames}_{start}_{end}"),
		REF = map2_chr(seqnames, start - 1, get_REF_entry),
		ALT = str_glue('<{CNV_type}>'),
	    QUAL = '.',
		FILTER = '.',
		INFO = str_glue("END={end};SVTYPE={CNV_type};SVLEN={length}"),
		# We can not phase the CNVs, so it's more accurate to set  ./.
		genotype = ifelse(CNV_type == 'LOH', '0/0', './.'),
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
				'##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">',
				'##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">')

vcf.format.content <- list(
	c('genotype', 'GT', "Segment genotype"),
	c('copynumber', 'CN', "Segment most-likely or estimated copy-number call"),
	c('numsnp', 'PN', "Number of points (i.e. SNP probes) in the segment"),
	c('CNV_caller', 'CA', "Segment called as CNV by these callers")#,
	#c('caller_confidence ', 'CC', "Caller dependent confidence score segment call (if available)")
)
vcf_cols <- c('CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT')

write_to_vcf <- function(tb, outvcf) {

	format_str <- paste(sapply(vcf.format.content, function (x) x[2]), collapse = ':')
	use.cols <- sapply(vcf.format.content, function (x) x[1])

	format.lines <- sapply(vcf.format.content, function (x) {q
		#Integer vs numeric ?!
		type <- ifelse(is.integer(tb[[x[1]]]), 'Integer', 'String')
		shorthand <- x[2]
		desc <- x[3]
		str_glue('##FORMAT=<ID={shorthand},Number=1,Type={type},Description="{desc}">')
	})

	writeLines(c(vcf.header, info.lines, format.lines,
	             paste0('#', paste(c(vcf_cols, sample_id), collapse = '\t'))), con = outvcf)

	tb <- tb %>%
		filter(CNV_type %in% args$include_state)

	if (nrow(tb) > 0) {
		tb %>%
			rowwise() %>%
			mutate(
				FORMAT = format_str,
			) %>%
			do({
				out = as.data.frame(.)
				out$sample_formatted = out[,use.cols] %>% paste(collapse = ':')
				out
			}) %>%
			dplyr::select(one_of(c(vcf_cols, 'sample_formatted', 'sample_id'))) %>%
			pivot_wider(names_from = sample_id, values_from = sample_formatted) %>%
			arrange(CHROM, POS) %>%
			write_tsv(outvcf, append=T)
	}

	}


if (args$mode == 'split-tools') {
	lapply(config$settings$CNV.calling.tools, function(use.tool) {
		outvcf <- str_glue('{data_path}/{sample_id}/{sample_id}.{use.tool}-cnv-calls{name_addition}.vcf')
		processed.calls %>%
			filter(caller_merging_state != 'combined' & CNV_caller == use.tool) %>%
			write_to_vcf(., outvcf = outvcf)
	})
} else {
	outvcf <- str_glue('{data_path}/{sample_id}/{sample_id}.combined-cnv-calls{name_addition}.vcf')
	processed.calls %>%
		filter(caller_merging_state != 'pre-overlap') %>%
		#some purrr walk function instead?
		write_to_vcf(., outvcf = outvcf)
}
