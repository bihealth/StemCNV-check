suppressMessages(library(tidyverse))
suppressMessages(library(optparse))
suppressMessages(library(DNAcopy))

parser <- OptionParser(
  usage = "usage: %prog [options] /path/to/inputfile.tsv /path/to/outputfile.tsv"
)

parser <- add_option(parser, c("-s", "--sd-undo"), 
                     default=1, help="Value for split SD undo")
 
args <- parse_args(parser, positional_arguments = 2)

sd.undo.val <- args$options$`sd-undo`
inputfile <- args$args[1]
outputfile <- args$args[2]

tb <- read_tsv(inputfile, show_col_types = FALSE) %>% 
  dplyr::select(-Index) %>%
  rename_with(~ str_remove(., '.*\\.'))

# tb$chip <- basename(args$args) %>% str_split('_') %>% .[1]
# tb$sample <- basename(args$args) %>% str_split('_') %>% .[2]
# tb$filter <- basename(args$args) %>% str_split('_') %>% .[3]

sampleID = basename(inputfile) %>% str_remove('\\..*$')

cna.basic <- CNA(tb$`Log R Ratio`, tb$Chr, tb$Position, data.type = 'logratio', sampleid = sampleID)
cna.basic.smoothed <- smooth.CNA(cna.basic)

cna.basic.smoothed.segmented <- segment(cna.basic.smoothed, min.width = 5, undo.splits = 'sdundo', undo.SD = sd.undo.val)

# p-values aren't much better - they are very sig. even for very short segments
#segments.p(cna.basic.smoothed.segmented)

segments.summary(cna.basic.smoothed.segmented) %>%
  write_tsv(outputfile)