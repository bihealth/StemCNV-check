# Readme

## About

StemCNV-check is a tool written to simplify copy number variation (CNV) analysis of SNP array data, specifically for quality control of (pluripotent) stem cell lines. 
StemCNV-check uses snakemake to run the complete analysis from raw data (.idat) up report generation for all defined samples with a single command. Samples need to be defined in a (tabular) sample table and the workflow settings are defined through a yaml file. 

StemCNV-check requires a linux environment (or WSL on windows) and a working conda (or mamba) installation. Follow the recommended instructions to install [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) or [mamba](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html).

## Installation

1. Clone the repository
2. Run `conda env create -f envs/base.yml` to create a base conda environment from which to run StemCNV-check
3. You may have to make the `StemCNV-check.py` script executable with `chmod +x StemCNV-check.py
4. All further dependencies (conda environments and docker containers) will be pulled automatically by snakemake when running the analysis

cnv-pipeline

## Setup

StemCNV-check requires a sample table and a config file to run. 

The sample table (default: sample_table.txt) is a tab-separated file describing all samples to be analyzed:
- Required columns: Sample_ID, Chip_Name, Chip_Pos, Sex, Reference_Sample
- Optional columns (reserved): Sample_Name, Regions_of_Interest
- See the `sample_table_example.txt` file for a description of individual columns
- Example files can be created using `StemCNV-check.py -a setup-files`

The config file (default: config.yaml) defines all settings for the analysis and inherits from the inbuilt default.  
Required settings that are not defined by default include static files specific to the used array platform and genome build:
- egt_cluster_file: the illumina cluster file (.egt) for the array platform, available from Illumina or the provider running the array 
- bpm_manifest_file: the beadpool manifest file (.bpm) for the array platform, available from Illumina or the provider running the array
- csv_manifest_file (optional): the manifest file in csv format, available from Illumina or the provider running the array

Additionally, the config file needs to define the following paths:
- raw_data_folder: path to the input directory under which the raw data (.idat) can be found. Ths folder should contain subfolders that match the Chip_Name column in the sample table (containing the array chip IDs)
- data_path: the output of StemCNV-check will be written to this path
- log_path: the log files of StemCNV-check will be written to this path


## Usage

Before the first analysis sample table and config file need to be set up (see above).
Automatic generation of the additional array & genome-build specific static files can only be done if sample data for 
that array is available.  
*Note*: unless provided directly this will also include download of fasta and gtf file for the reference genome build.

`StemCNV-check.py -a make-staticdata --genome <hg19|hg38> --snp-array-name <array_platform> [-s <sample_table> -c <config_file>]`

To run the analysis:

`StemCNV-check.py [-s <sample_table> -c <config_file> -a run]`

## Example data

This repository contains example data (using data from the Genome in a Bottle samples) that can be used to test the setup.
After pulling the repository and creating and activating the base scnc-quac conda environment, run the following commands 
(Note that this will download a fasta and gtf file for the human genome. If you have suitable files available locally, 
it is recommended to replace the corresponding paths in the config.yaml to avoid unnecessary and time-consuming downloads):
- `cd test_data`
- `../StemCNV-check.py -a make-staticdata --genome hg19 --snp-array-name testrun` 
- `../StemCNV-check.py`

## Output

StemCNV-check will produce the following output files for each sample, when run with default settings:
- `data_path/{sample}/{sample}.unprocessed.vcf`; `data_path/{sample}/{sample}.processed-data.tsv`  
  The unfiltered SNP data of the array in vcf and tabular format (vcf contains more information)
- `data_path/{sample}/{sample}.stats.txt`  
  The CNV calls for the sample GenCall stat
- `data_path/{sample}/{sample}.filtered-data-{filter}.tsv`  
  Filtered SNP data in tabular format, the default filter is 'extended' (see default config), other filters can be defined
- `data_path/{sample}/{sample}.CBS.tsv`  
  The CNV calls for the sample from the CBS (Circular Binary Segmentation) algorithm
- `data_path/{sample}/{sample}.penncnv-{auto|chrx|chry}.tsv`  
  The CNV calls for the sample from the PennCNV caller (Note: the chry file will only be created for samples annotated as male)
- `data_path/{sample}/{sample}.combined-cnv-calls.{tsv|vcf}`  
  The CNV calls processed, combined and annotated by StemCNV-check in tabular and vcf format. Annotation includes comparison against reference sample, call scoring and gene annotation.
- `data_path/{sample}/{sample}.user-report.html`; ...  
  Html reports containing summary statistics, QC statistics, lists CNV calls sorted by annotation score, plots of most/all CNVs and sample comparison. The default 'user-report' only contains plots for critical and reportable calls, a full report can easily be enabled in the config.yaml. The content of either the default or any additional reports can also be fine-tuned through the config.yaml file.
- `data_path/{sample}/{sample}.summary-check.tsv`  
  A tabular export of the summary table contained in the report