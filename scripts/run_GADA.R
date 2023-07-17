#! /usr/bin/Rscript
# Run GADA & MAD CNV calling
suppressMessages(library(argparse))

parser <- ArgumentParser(description="Run GADA and MAD CNV calling")

parser$add_argument('inputfile', type = 'character', help='Path to input file')
parser$add_argument('outputfile', type = 'character', help='Path to output file')

parser$add_argument("-a", "--aalpha", type = 'numeric', default = 0.2, help="aAlpha for SBL function")
parser$add_argument("-g", "--gadat", default = 15, help="T for BE of GADA")
parser$add_argument("-m", "--madt", default = 5, help="T for BE of MAD")
parser$add_argument("-l", "--minseglen", default = 5, help="Min Probe number for calls")

args <- parser$parse_args()
#args <- parser$parse_args(c('/home/vonkunic_c/Misc-Projects/CNV-pipeline/test/data/BIHi005-A-13/BIHi005-A-13.filtered-data.full.tsv', '/home/vonkunic_c/Misc-Projects/CNV-pipeline/test/data/BIHi005-A-13/BIHi005-A-13.GADA.full.tsv'))

suppressMessages(library(tidyverse))
suppressMessages(library(gada))
suppressMessages(library(mad))
suppressMessages(library(plyranges))


# TODO might want to keep / capture the extra cols & printed messages for tool QC stats ?

inputfile <- args$inputfile
outputfile <- args$outputfile

aalpha <- args$aalpha
gada.T <- args$gadat
mad.T <- args$madt
minseglen <- args$minseglen

#GADA will have issues if some expexted chromosomes don't have data
#Maybe this can be pulled from config?
tb <- read_tsv(inputfile, show_col_types = FALSE)
all.chr <- c(1:22, 'X', 'Y')
use.chr <- all.chr[all.chr %in% unique(tb$Chr)]

sampleID <- basename(inputfile) %>% str_remove('\\.filtered-data\\..*\\.tsv$')

gada.obj <- setupGADA(inputfile,
							MarkerIdCol=2,
							ChrNameCol=4, ChrPosCol=5,
							log2ratioCol=16, BAFcol=15,
							chrs = use.chr)

gen.info <- attr(gada.obj, 'gen.info')
# --------------------------------------------------------------------
# 	(higher sensitivity , higher FDR ) <----> $(\alpha = 0.2,T > 3)$
# 									   <----> $(\alpha = 0.5,T > 4)$
# 	(lower sensitivity , lower FDR )   <----> $(\alpha = 0.8,T > 5)$
# -------------------------------------------------------------------

gada.obj <- SBL(gada.obj, estim.sigma2 = TRUE, aAlpha = aalpha)
gada.obj <- BackwardElimination(gada.obj, T=gada.T, MinSegLen=minseglen)

cnvs <- summary(gada.obj, print=FALSE) %>%
	# this pots some useful info to log ?
	as_tibble() %>%
	filter(State != 0) %>%
	mutate(State = ifelse(State == 1, 'gain', 'loss'),
		   chromosome = paste0('chr', chromosome)) %>%
	as_granges(seqnames = chromosome,
				 start = IniProbe, end = EndProbe,
				 keep_mcols = T) %>%
	select(-MeanAmp)

# Only the 'par' MAD functions seem to work properly
tmp <- tempdir()
oldwd <- getwd()
setwd(tmp)
dir.create("rawData")
file.symlink(inputfile, file.path(tmp, 'rawData', basename(inputfile)))

if (!all(all.chr %in% use.chr)) {
	warning('MAD can only calculated mosaicism for both sex chromomsomes together, at least one is missing (Normal for filtered female data).')
	all.chr <- 1:22
	use.chr <- all.chr[all.chr %in% use.chr]
}

mad.col <- setupParGADA.B.deviation(NumCols=ncol(tb), GenoCol=11,
									BAFcol=15, log2ratioCol=16,
									MarkerIdCol=2,
									ChrNameCol=4, ChrPosCol=5,
									chrs = use.chr)

parSBL(mad.col, estim.sigma2=TRUE, aAlpha=aalpha)

parBE.B.deviation(mad.col, T=mad.T, MinSegLen=minseglen)

# The number codes correspond to the following abnormalities: UPD (1), deletion (2), duplication (3), trisomy (4) and LOH (5)
# 1 UPD = uniparental disomy == LOH
mosaic.states <- c('UPD', 'deletion', 'duplication', 'trisomy', 'LOH')
mosaics <- tryCatch(getMosaics(mad.col) %>%
						mutate(State = mosaic.states[State]) %>%
						select(-LRR, -LRR.se, -Bdev, -sample),
				   error = function(e) {
					   if (str_detect(as.character(e), "stop_if_wrong_length")) {
						   message('No moasics detected with MAD. Try decreasing mad.T value.')
					   } else {
						   stop(e)
					   }
					   GRanges() }
)

setwd(oldwd)

calls.out <- bind_ranges(cnvs, mosaics) %>%
	as_tibble() %>%
	dplyr::rename(numsnp = LenProbe, CNV.state = State,
				  length = width, Chr = seqnames) %>%
	mutate(sample_id = sampleID,
	       tool = 'GADA') %>%
	dplyr::select(-strand)

write_tsv(calls.out, outputfile)
