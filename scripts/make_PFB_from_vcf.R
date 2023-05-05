suppressMessages(library(tidyverse))
suppressMessages(library(optparse))

# Script to make PFB file based on information in vcf file (derived from Illumina cluster file)

parser <- OptionParser(
	usage = "usage: %prog /path/to/inputfile.vcf /path/to/outputfile.pfb"
)

args <- parse_args(parser, positional_arguments = 2)

inputfile <- args$args[1]
outputfile <- args$args[2]

#TODO: use vcfR instead

vcf.info <- read_delim(inputfile, comment = '#', delim = '\t', 
                              col_names = c('Chr', 'Position', 'Name', 'Ref', 'Alt', 'x1', 'x2', 'INFO', 'FORMAT', 'SAMPLE')) %>%
  select(Name, Chr, Position, Ref, Alt, INFO) %>%
  separate(INFO, 
           .$INFO[[1]] %>% str_remove_all('=[0-9.]+') %>% str_split(';') %>% unlist(),
           sep=';[^=]+=', convert = T) %>%
  mutate(Chr = str_remove(Chr, 'chr'),
         GC = str_remove(GC, 'GC='),
         alleles = paste0(Ref, ',', Alt),
  			 A_freq = (2*N_AA + N_AB) / (2*(N_AA + N_AB + N_BB)),
  			 #This is B_allele frequency, or 'population frequency B-allele'. PennCNV needs this
  			 PFB = (2*N_BB + N_AB) / (2*(N_AA + N_AB + N_BB)),
  			 #This is AF/Alternate allele frequency. Some bcftool modules need this
  			 Alt_freq = ifelse(ALLELE_B == 1, PFB, A_freq)
  			 ) 

vcf.info %>%
  select(Name, Chr, Position, PFB) %>%
	mutate(Chr = factor(Chr, levels=c(1:22, 'X', 'Y'))) %>%
  filter(if_all(.fns = ~!is.na(.))) %>%
	arrange(Chr, Position) %>%
  write_tsv(outputfile)
