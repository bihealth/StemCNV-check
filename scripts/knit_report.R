suppressMessages(library(knitr))
suppressMessages(library(tidyverse))
suppressMessages(library(optparse))

parser <- OptionParser(
  usage = "usage: %prog [options] <sample_id> <basepath> <snakedir> /path/to/sample_table.tsv /path/to/config.yaml"
)

# parser <- add_option(parser, c("-n", "--no-reference"), action="store_true",
#                      default=FALSE, help="Ignore reference sample regardless of samplesheet")

args <- parse_args(parser, positional_arguments = 5)


#TODO make all paths absolute

sample_id   <- args$args[1]
basepath    <- args$args[2] %>% normalizePath
snakedir    <- args$args[3] %>% normalizePath
sampletable <- args$args[4] %>% normalizePath
configfile  <- args$args[5] %>% normalizePath

reportfile  <- file.path(snakedir, "scripts", "CNV_report.Rmd") %>% normalizePath
workdir     <- file.path(basepath, 'data', sample_id) %>% normalizePath
outfile     <- paste0(sample_id, ".CNV-report.html") #%>% normalizePath

# clear previously generated images
if (dir.exists(file.path(workdir, 'report_images'))) {
	system(str_glue('rm {workdir}/report_images/*'))
}

# Run Rmd - all files should be stored in the final output folder (== workdir)
setwd(workdir)
rmarkdown::render(
  input = reportfile, 
  output_dir = workdir,
  intermediates_dir = workdir,
  knit_root_dir = workdir,
  output_file = outfile,
  output_options = list(self_contained = TRUE), # base.dir = figure.dir ?
  params = list(base_path = basepath,
                sample_id = sample_id,
                #reference_id = reference_id,
                sampletable = sampletable,
                configfile = configfile,
                workdir = workdir
                )
)
