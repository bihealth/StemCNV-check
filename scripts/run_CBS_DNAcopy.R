#! /usr/bin/Rscript
# Run LRR segmentation with CBS/DNACopy
suppressMessages(library(argparse))

parser <- ArgumentParser(description="Run LRR segmentation with CBS/DNACopy")

parser$add_argument('inputfile', type = 'character', help='Path to input file')
parser$add_argument('outputfile', type = 'character', help='Path to output file')
parser$add_argument('configfile', type = 'character', help='Path to config file')
parser$add_argument('sexfile', type = 'character', help='Path to sexfile')

parser$add_argument('-s', '--sd-undo', type = 'numeric', default = 1,
					help="Value for split SD undo")

args <- parser$parse_args()
#args <- parser$parse_args(c('/home/vonkunic_c/Misc-Projects/CNV-pipeline/test/data/BIHi005-A/BIHi005-A.filtered-data.full.tsv', 'test.out', 'default_config.yaml', 'penncnv-sexfile.txt'))

suppressMessages(library(tidyverse))
suppressMessages(library(DNAcopy))
suppressMessages(library(yaml))

sd.undo.val <- args$sd_undo
inputfile <- args$inputfile
outputfile <- args$outputfile
config <- read_yaml(args$configfile)

sextable <- read_tsv(args$sexfile, col_names=c('filename', 'sex'))
sex <- sextable[sextable$filename == inputfile, ]$sex

# CBS gain/loss
CBS.LRR.th.value 	  <- config$settings$CBS$LRR.th.value
CBS.LRR.th.value.Xadj <- config$settings$CBS$LRR.th.value.Xadj
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

sampleID <- basename(inputfile) %>% str_remove('\\.filtered-data\\..*\\.tsv$')

cna.basic <- CNA(tb$`Log R Ratio`, tb$Chr, tb$Position, data.type = 'logratio', sampleid = sampleID)
cna.basic.smoothed <- smooth.CNA(cna.basic)

#TODO: parameterise min.width
cna.basic.smoothed.segmented <- segment(cna.basic.smoothed, min.width = 5, undo.splits = 'sdundo', undo.SD = sd.undo.val)

# p-values aren't much better - they are very sig. even for very short segments
#segments.p(cna.basic.smoothed.segmented)

tb <- segments.summary(cna.basic.smoothed.segmented) %>%
  		dplyr::rename(Chr = chrom, start = loc.start, end = loc.end,
									numsnp = num.mark, sample_id = ID) %>%
	rowwise %>%
	mutate(
		sample_id = sampleID,
		length = end - start, # TODO: open / half open / +- 1 ??
		Chr = paste0('chr', Chr),
		Chr = factor(Chr, levels = c(paste0('chr', 1:22), 'chrX', 'chrY')),
		snp.density = numsnp / length * 1e6,
		copynumber = case_when(
			#Male is default CN=1 on X & Y, also unqiue cutoffs
			sex == 'm' & Chr %in% c('chrX', 'chrY') & seg.median < LRR.male.XorY.loss       ~ 0,
			sex == 'm' & Chr %in% c('chrX', 'chrY') & seg.median > LRR.male.XorY.gain.large ~ 3,
			sex == 'm' & Chr %in% c('chrX', 'chrY') & seg.median > LRR.male.XorY.gain       ~ 2,
			# Unique cutoffs for female X (behaves different from autosome)
			sex == 'f' & Chr == 'chrX' & seg.median < LRR.female.XX.loss                    ~ 0,
			sex == 'f' & Chr == 'chrX' & seg.median < LRR.female.X.loss                     ~ 1,
			sex == 'f' & Chr == 'chrX' & seg.median > LRR.female.X.gain.large               ~ 4,
			sex == 'f' & Chr == 'chrX' & seg.median > LRR.female.X.gain                     ~ 3,
			# Default cutoffs
			seg.median < LRR.loss.large                                                     ~ 0,
			seg.median < LRR.loss                                                           ~ 1,
			seg.median > LRR.gain.large                                                     ~ 4,
			seg.median > LRR.gain                                                           ~ 3,
			TRUE ~ 2
		),
		CNV.state = case_when(# Male is default CN=1 on X & Y
			sex == 'm' & Chr %in% c('chrX', 'chrY') & copynumber == 2     ~ 'gain',
			sex == 'm' & Chr %in% c('chrX', 'chrY') & copynumber == 1     ~ NA,
			copynumber > 2                                                ~ 'gain',
			copynumber < 2                                                ~ 'loss ',
			.default = NA
		),
		tool = 'CBS',
		conf = NA,
		# Older version
		CNV.state.old = ifelse(seg.median < -CBS.LRR.th.value, 'loss', NA),
		CNV.state.old = ifelse(seg.median > CBS.LRR.th.value, 'gain', CNV.state.old),
		CNV.state.old = ifelse(Chr == 'chrX', NA, CNV.state.old),
		CNV.state.old = ifelse(Chr == 'chrX' & seg.median < -CBS.LRR.th.value + ifelse(sex == 'm', -1, 1) * CBS.LRR.th.value.Xadj, 'loss', CNV.state.old),
		CNV.state.old = ifelse(Chr == 'chrX' & seg.median > CBS.LRR.th.value + ifelse(sex == 'm', -1, 1) * CBS.LRR.th.value.Xadj, 'gain', CNV.state.old),
		copynumber.old = ifelse(is.na(CNV.state.old), 2, 3),
		copynumber.old = ifelse(CNV.state.old == 'loss', 1, copynumber.old),
	) %>%
	filter(!is.na(CNV.state) & !is.na(Chr))

if (sex == 'f') {
	tb <- filter(tb, Chr != 'chrY')
}

write_tsv(tb, outputfile)
