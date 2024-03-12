# Readme

## About

<...> is a tool written to simplify copy number variation (CNV) analysis of SNP array data, specifically for quality control of (pluripotent) stem cell lines. 
<...> uses snakemake to run the complete analysis from raw data (.idat) up report generation for each sample with a single command.

It requires a linux environment (including WSL) and a working conda (or mamba) installation.

## Installation

1. Clone the repository
2. Run `conda env create -f envs/base.yml` to create a base conda environment from which to run <...>
3. All further dependencies (conda environments and docker containers) will be pulled automatically by snakemake when running the analysis

## Setup

<...> requires a sample table and a config file to run. 
sample_table.txt (tab-separated)

- sample_id
- sex
- ref
- ...

The config file (default: config.yaml) needs to define static files specific to the used array platform and genome build:
- egt
- bpm
- csv
- gtf
- genome.fa

Additionally, the config file needs to define the following paths:
- input_dir: path to the directory containing the raw data (.idat files)
- output_dir: path to the directory where the analysis results will be stored
- log_dir: path to the directory where the log files will be stored

Examples of these files can be created using this command: `<...> -a setup-files`

## Usage

Before the first analysis additional array & genome-build specific files need to be created:

`<...> -a make-staticdata --genome <genome_build> --snp-array-name <array_platform> [-s <sample_table> -c <config_file>]`

To run the analysis:

`<...> [-s <sample_table> -c <config_file> -a run]`

## Output

