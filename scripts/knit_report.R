suppressMessages(library(knitr))
suppressMessages(library(tidyverse))
suppressMessages(library(optparse))

parser <- OptionParser(
  usage = "usage: %prog [options] <sample_id> <basepath> <snakedir> /path/to/sample_table.tsv /path/to/config.yaml"
)

parser <- add_option(parser, c("-n", "--no-reference"), action="store_true",
                     default=FALSE, help="Ignore reference sample regardless of samplesheet")

args <- parse_args(parser, positional_arguments = 5)

sample_id   <- args$args[1]
basepath    <- args$args[2]
snakedir    <- args$args[3]
sampletable <- args$args[4]
configfile  <- args$args[5]

reportfile  <- file.path(snakedir, "scripts", "CNV_report.Rmd")
workdir     <- file.path(basepath, 'data', sample_id)
outfile     <- paste0(sample_id, ".CNV-report.html")

samples.tb <- read_tsv(sampletable, col_types = 'cccccc')
reference <- samples.tb[samples.tb$Sample_ID == sample_id, ]$ReferenceSample
reference_id <- ifelse(is.na(reference), '', samples.tb[samples.tb$SampleName == reference, ]$Sample_ID)

setwd(workdir)

# message(getwd())
# message(reference_id)

rmarkdown::render(
  input = reportfile, 
  output_dir = workdir,
  intermediates_dir = workdir,
  knit_root_dir = workdir,
  output_file = outfile,
  output_options = list(self_contained = TRUE), # base.dir = figure.dir ?
  params = list(base_path = basepath,
                sample_id = sample_id,
                reference_id = reference_id,
                sampletable = sampletable,
                configfile = configfile,
                workdir = workdir
                )
)
