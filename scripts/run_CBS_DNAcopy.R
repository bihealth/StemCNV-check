suppressMessages(library(tidyverse))
suppressMessages(library(optparse))
suppressMessages(library(DNAcopy))

parser <- OptionParser(
  usage = "usage: %prog [options] /path/to/inputfile.tsv /path/to/outputfile.tsv"
)

parser <- add_option(parser, c("-s", "--sd-undo"), 
                     default=1, help="Value for split SD undo")
 
args <- parse_args(parser, positional_arguments = 2)
#args <- parse_args(parser, positional_arguments = 2, c('data/205857150048_R01C01/205857150048_R01C01.filtered-data.full.tsv', 'test.out'))

sd.undo.val <- args$options$`sd-undo`
inputfile <- args$args[1]
outputfile <- args$args[2]

tb <- read_tsv(inputfile, show_col_types = FALSE) %>% 
  dplyr::select(-Index) %>%
  rename_with(~ str_remove(., '.*\\.'))

# tb$chip <- basename(args$args) %>% str_split('_') %>% .[1]
# tb$sample <- basename(args$args) %>% str_split('_') %>% .[2]
# tb$filter <- basename(args$args) %>% str_split('_') %>% .[3]

#TODO: XY probes -> separate or remove?

#TODO:
# add tresholds etc here
# # CBS gain/loss
# CBS.LRR.th.value 	  <- config$settings$CBS$LRR.th.value
# CBS.LRR.th.value.Xadj <- config$settings$CBS$LRR.th.value.Xadj

sampleID <- basename(inputfile) %>% str_remove('\\..*$')

cna.basic <- CNA(tb$`Log R Ratio`, tb$Chr, tb$Position, data.type = 'logratio', sampleid = sampleID)
cna.basic.smoothed <- smooth.CNA(cna.basic)

cna.basic.smoothed.segmented <- segment(cna.basic.smoothed, min.width = 5, undo.splits = 'sdundo', undo.SD = sd.undo.val)

# p-values aren't much better - they are very sig. even for very short segments
#segments.p(cna.basic.smoothed.segmented)

segments.summary(cna.basic.smoothed.segmented) %>%
  	# 	dplyr::rename(Chr = chrom, start = loc.start, end = loc.end,
	# 								numsnp = num.mark, sample_id = ID) %>%
	# 	# TODO -> most of this should probably be moved into the CBS script
	# 	mutate(
	# 		sample_id = str_remove(sample_id, '^X'),
	# 		length = end - start, # TODO: open / half open / +- 1 ??
	# 		Chr = paste0('chr', Chr),
	# 		Chr = factor(Chr, levels = c(paste0('chr', 1:22), 'chrX', 'chrY')),
	# 		snp.density = numsnp / length * 1e6,
	# 		CNV.state = ifelse(seg.median < -CBS.LRR.th.value, 'loss', NA),
	# 		CNV.state = ifelse(seg.median > CBS.LRR.th.value, 'gain', CNV.state),
	# 		CNV.state = ifelse(Chr == 'chrX', NA, CNV.state),
	# 		CNV.state = ifelse(Chr == 'chrX' & seg.median < -CBS.LRR.th.value + ifelse(sex == 'm', -1, 1) * CBS.LRR.th.value.Xadj, 'loss', CNV.state),
	# 		CNV.state = ifelse(Chr == 'chrX' & seg.median > CBS.LRR.th.value + ifelse(sex == 'm', -1, 1) * CBS.LRR.th.value.Xadj, 'gain', CNV.state),
	# 		tool = 'CBS',
	# 		conf = NA,
	# 		#TODO could try double / 3x / ... cutoff for 0/4 CN ?
	# 		# > 0 doesn't make sense though (would mean no detection == filtered out)
	# 		# > 4 maybe? very hard to tell & probably irrelevant
	# 		copynumber = ifelse(is.na(CNV.state), 2, 3),
	# 		copynumber = ifelse(CNV.state == 'loss', 1, copynumber),
	# 	) %>%
	# 	filter(!is.na(CNV.state) & !is.na(Chr))
	# if (sex == 'f') {
	# 	tb <- filter(tb, Chr != 'chrY')
	# }
  write_tsv(outputfile)