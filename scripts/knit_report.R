#! /usr/bin/Rscript
# Wrapper to knit Rmd report
suppressMessages(library(argparse))

parser <- ArgumentParser(description="Wrapper to knit Rmd report")

parser$add_argument('sample_id', type = 'character', help='Sample_ID to use')
parser$add_argument('configfile', type = 'character', help='Path to config file')

args <- parser$parse_args()

suppressMessages(library(knitr))
suppressMessages(library(tidyverse))
suppressMessages(library(yaml))

sample_id   <- args$sample_id
configfile  <- args$configfile %>% normalizePath

config <- read_yaml(configfile)

basepath    <- config$basedir %>% normalizePath
snakedir    <- config$snakedir %>% normalizePath
sampletable <- config$sample_table %>% normalizePath
datapath    <- config$data_path

reportfile  <- file.path(snakedir, "scripts", "CNV_report.Rmd") %>% normalizePath
workdir     <- file.path(basepath, datapath, sample_id) %>% normalizePath
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
