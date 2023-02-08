How to create stattic files:

# PFB

## Option 1) PennCNV - compile from existing data
-> requires a relatively large set of samples, ideally with high heterogeneity
First run pipeline to get processed data from all samples you have access to then use this:

`$CONDA_PREFIX/pipeline/PennCNV-1.0.5/compile_pfb.pl data/*/*.processed-data.tsv -out static-data/pfb_compiled_from_PennCNV.pfb`

## Option 2) Extra allele frequencies from the info of the egt_cluster file, which is written into vcf's (parsing egt file is harder)
First run pipeline on any given sample to get a vcf file

TODO: make the script cmd line callable (extract from the other one)

`Rscript make_PFB_from_VCF.R <vcf-file> static-data/pfb_extracted_from_egt.pfb`


# Allele Frequency files (not needed anymore?)
...


# GCmodel files

## Option 1) PennCNV script
PennCNV ships with UCSC GC annotation files, that contain GC content in (?) ~5kb windowes
Using the `cal_gc_snp.p` script these files can be used to calculate GC content in 1MB windows surrounding each marker

`$CONDA_PREFIX/pipeline/PennCNV-1.0.5/cal_gc_snp.pl $CONDA_PREFIX/pipeline/PennCNV-1.0.5/gc_file/hg38.gc5Base.txt static-data/pfb_extracted_from_egt.pfb -out static-data/gcmodel_hg38PennCNV_pfbEGTaf.tsv`

## Option 2) Write your own script to calculate GC from genome fasta
TODO: use some suitable R or python library ...