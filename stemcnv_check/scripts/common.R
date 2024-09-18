if (snakemake@config$verbose_level >= 2) {
    save.image(file = file.path(dirname(snakemake@log[['err']]), 'R_env.Rdata'))
}
log <- file(snakemake@log[['err']], 'wt')
sink(log, append = T)
sink(log, append = T, type = 'message')
