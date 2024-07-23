#! /usr/bin/Rscript
# Wrapper to knit Rmd report
suppressMessages(library(argparse))

parser <- ArgumentParser(description="Wrapper to knit Rmd report")

parser$add_argument('sample_id', type = 'character', help='Sample_ID to use')
parser$add_argument('report_name', type = 'character', help='Which report (defined in config) to run')
parser$add_argument('configfile', type = 'character', help='Path to config file')

args <- parser$parse_args()

suppressMessages(library(knitr))
suppressMessages(library(tidyverse))
suppressMessages(library(yaml))

sample_id   <- args$sample_id
report_name <- args$report_name
configfile  <- args$configfile %>% normalizePath

config <- read_yaml(configfile)

if (!report_name %in% names(config$reports)) {
  stop(str_glue('Error: {report_name} ir not a defined report in the config: {configfile}'))
}
report_settings <- config$reports[[report_name]]
filetype    <- report_settings$file_type

basepath    <- config$basedir %>% normalizePath
snakedir    <- config$snakedir %>% normalizePath
sampletable <- config$sample_table %>% normalizePath
datapath    <- config$data_path

report.template  <- file.path(snakedir, "scripts", "report_template.Rmd") %>% normalizePath
outfile     <- str_glue("{sample_id}.{report_name}.{filetype}")

version <- readLines(file.path(snakedir, 'version.py')) %>% str_trim() %>%
  str_remove('.* = ') %>% str_remove_all("\'")

if (fs::is_absolute_path(datapath)) {
  workdir   <- file.path(datapath, sample_id) %>% normalizePath
} else {
  workdir   <- file.path(basepath, datapath, sample_id) %>% normalizePath
}
# clear previously generated images
if (dir.exists(str_glue('{workdir}/{report_name}_images'))) {
	system(str_glue('rm {workdir}/{report_name}_images/*'))
}

# Run Rmd - all files should be stored in the final output folder (= workdir)
# Note: the intermediates_dir/knit_root_dir settings only apply do (r)md, but *not* to latex for pdf generation (https://github.com/rstudio/rmarkdown/issues/1975)
# In order to have latex output everything (incl error logs) in the 'output_dir' the input template would need to (temporarily) copied over
setwd(workdir)
rmarkdown::render(
  input = report.template,
  output_dir = workdir,
  intermediates_dir = workdir,
  knit_root_dir = workdir,
  output_format = paste0(filetype, '_document'),
  output_file = outfile,
  output_options = list(self_contained = TRUE),
  params = list(base_path = basepath,
                sample_id = sample_id,
                sampletable = sampletable,
                report_name = report_name,
                out_format = filetype,
                configfile = configfile,
                workdir = workdir,
                version = version
                )
)
