#! /usr/bin/Rscript
# Run LRR segmentation with CBS/DNACopy
suppressMessages(library(argparse))

parser <- ArgumentParser(description="Run LRR segmentation with CBS/DNACopy")

parser$add_argument('inputfile', type = 'character', help='Path to input file')
parser$add_argument('outputfile', type = 'character', help='Path to output file')
parser$add_argument('configfile', type = 'character', help='Path to config file')
parser$add_argument('sampletable', type = 'character', help='Path to sampletable')

parser$add_argument('-s', '--sd-undo', type = 'double', default = 1,
					help="Value for split SD undo")

args <- parser$parse_args()

suppressMessages(library(tidyverse))
suppressMessages(library(DNAcopy))
suppressMessages(library(yaml))

inputfile <- args$inputfile
outputfile <- args$outputfile
config <- read_yaml(args$configfile)
sd.undo.val <- args$sd_undo
min.width <- config$settings$CBS$min.width

sampleID <- basename(inputfile) %>% str_remove('\\.filtered-data-.*\\.tsv$')

sampletable <- read_tsv(args$sampletable, col_types = 'cccccc', comment = '#')
sex <- sampletable[sampletable$Sample_ID == sampleID, ]$Sex %>%
	tolower() %>% substr(1, 1)

# More fine grained loss/gain thresholds, based on Samules Nx settings
LRR.loss <- config$settings$CBS$LRR.loss
LRR.loss.large <- config$settings$CBS$LRR.loss.large
LRR.gain <- config$settings$CBS$LRR.gain
LRR.gain.large <- config$settings$CBS$LRR.gain.large
# plus specifics for sex chromosomes
LRR.male.XorY.loss <- config$settings$CBS$LRR.male.XorY.loss
LRR.male.XorY.gain <- config$settings$CBS$LRR.male.XorY.gain
LRR.male.XorY.gain.large <- config$settings$CBS$LRR.male.XorY.gain.large
LRR.female.X.loss <- config$settings$CBS$LRR.female.X.loss
LRR.female.XX.loss <- config$settings$CBS$LRR.female.XX.loss
LRR.female.X.gain <- config$settings$CBS$LRR.female.X.gain
LRR.female.X.gain.large <- config$settings$CBS$LRR.female.X.gain.large

tb <- read_tsv(inputfile, show_col_types = FALSE) %>% 
  dplyr::select(-Index) %>%
  rename_with(~ str_remove(., '.*\\.'))

cna.basic <- CNA(tb$`Log R Ratio`, tb$Chr, tb$Position, data.type = 'logratio', sampleid = sampleID)
cna.basic.smoothed <- smooth.CNA(cna.basic)

cna.basic.smoothed.segmented <- segment(cna.basic.smoothed, min.width = min.width, undo.splits = 'sdundo', undo.SD = sd.undo.val)

tb <- segments.summary(cna.basic.smoothed.segmented) %>%
  		dplyr::rename(Chr = chrom, start = loc.start, end = loc.end,
									numsnp = num.mark, sample_id = ID) %>%
	mutate(
		sample_id = sampleID,
		length = end - start,
		Chr = paste0('chr', Chr),
		Chr = factor(Chr, levels = c(paste0('chr', 1:22), 'chrX', 'chrY')),
		snp.density = numsnp / length * 1e6,
		copynumber = case_when(
			#Male is default CN=1 on X & Y, also unqiue cutoffs
			sex == 'm' & Chr %in% c('chrX', 'chrY') & seg.median < LRR.male.XorY.loss       ~ 0,
			sex == 'm' & Chr %in% c('chrX', 'chrY') & seg.median > LRR.male.XorY.gain.large ~ 3,
			sex == 'm' & Chr %in% c('chrX', 'chrY') & seg.median > LRR.male.XorY.gain       ~ 2,
			sex == 'm' & Chr %in% c('chrX', 'chrY')                                         ~ 1,
			# Unique cutoffs for female X (behaves different from autosome)
			sex == 'f' & Chr == 'chrX' & seg.median < LRR.female.XX.loss                    ~ 0,
			sex == 'f' & Chr == 'chrX' & seg.median < LRR.female.X.loss                     ~ 1,
			sex == 'f' & Chr == 'chrX' & seg.median > LRR.female.X.gain.large               ~ 4,
			sex == 'f' & Chr == 'chrX' & seg.median > LRR.female.X.gain                     ~ 3,
			sex == 'f' & Chr == 'chrX'                                                      ~ 2,
			# Default cutoffs
			seg.median < LRR.loss.large                                                     ~ 0,
			seg.median < LRR.loss                                                           ~ 1,
			seg.median > LRR.gain.large                                                     ~ 4,
			seg.median > LRR.gain                                                           ~ 3,
			TRUE ~ 2,
			.default = 2
		),
		CNV.state = case_when(# Male is default CN=1 on X & Y
			sex == 'm' & Chr %in% c('chrX', 'chrY') & copynumber == 2     ~ 'gain',
			sex == 'm' & Chr %in% c('chrX', 'chrY') & copynumber == 1     ~ NA,
			copynumber > 2                                                ~ 'gain',
			copynumber < 2                                                ~ 'loss ',
			.default = NA
		),
		tool = 'CBS',
		ID = paste(tool, CNV.state, Chr, start, end, sep='_'),
		tool_confidence = NA,
	) %>%
	filter(!is.na(CNV.state) & !is.na(Chr))

if (sex == 'f') {
	tb <- filter(tb, Chr != 'chrY')
}

write_tsv(tb, outputfile)
