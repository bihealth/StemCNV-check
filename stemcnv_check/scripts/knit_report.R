# Redirect warnings & errors to snakemake logging, save R environment if debugging
source(file.path(snakemake@config$snakedir, 'scripts/common.R'))

library(tidyverse)
library(knitr)

config <- snakemake@config
report_name <- snakemake@wildcards$report

if (!report_name %in% names(config$reports)) {
    stop(str_glue('Error: {report_name} is not a defined report in the config'))
}

report.template  <- file.path(config$snakedir, "scripts", "report_template.Rmd") %>% normalizePath
outfile     <- snakemake@output$report
rmd_workdir <- dirname(outfile) %>% normalizePath

# # clear previously generated images
# if (dir.exists(str_glue('{rmd_workdir}/{report_name}_images'))) {
# 	system(str_glue('rm {rmd_workdir}/{report_name}_images/*'))
# }

# Run Rmd - all files should be stored in the final output folder (= workdir)
# Note: the intermediates_dir/knit_root_dir settings only apply do (r)md, but *not* to latex for pdf generation (https://github.com/rstudio/rmarkdown/issues/1975)
# In order to have latex output everything (incl error logs) in the 'output_dir' the input template would need to (temporarily) copied over
# setwd(rmd_workdir)
rmarkdown::render(
    input = report.template,
    output_dir = rmd_workdir,
    intermediates_dir = rmd_workdir,
    knit_root_dir = rmd_workdir,
    output_format = paste0(snakemake@wildcards$ext, '_document'),
    output_file = outfile,
    output_options = list(self_contained = TRUE),
    params = list(
        sample_id =  snakemake@wildcards$sample_id,
        inputs = snakemake@input,
        config = config,
        report_config = snakemake@params$report_config,
        config_delta = snakemake@params$config_delta,
        #report_name = report_name,
        out_format = snakemake@wildcards$ext,
        #workdir = rmd_workdir,
        basedir = config$basedir,
        version = snakemake@params$version
    )
)
