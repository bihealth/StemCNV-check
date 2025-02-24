# Readme

## About

StemCNV-check is a tool written to simplify copy number variation (CNV) analysis of SNP array data, specifically for quality control of (pluripotent) stem cell lines. 
StemCNV-check uses snakemake to run the complete analysis from raw data (.idat) up report generation for all defined samples with a single command. Samples need to be defined in a (tabular) sample table and the workflow settings are defined through a yaml file. 

StemCNV-check requires a linux environment (or WSL on windows) and a working conda (or mamba) installation. Follow the recommended instructions to install [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) or [mamba](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html).

## Installation

StemCNV-check will be made available through the bioconda channel in the future. 

For now, only installation 'from source' is possible:

1. Clone this git repository
2. *optional, but recommended* Create a new enviroment, i.e. conda create -n stemcnv-check python=3.12, then activate it
   - Note: on some systems like WSL you may also need: `apptainer` and `gcc_linux-64` (<14 for recent datrie issue)
3. Install dependencies and the stemcnv-check script using pip `pip install -e .`. For development, use `pip install -e .[all]`
4. All further runtime dependencies (conda environments and docker containers) will be pulled automatically by snakemake when running the analysis

## Setup

StemCNV-check requires a sample table and a config file to run. Example files can be created using `stemcnv-check setup-files`.

The sample table (default: sample_table.tsv) is a tab-separated file describing all samples to be analyzed:
- Required columns: Sample_ID, Chip_Name, Chip_Pos, Array_Name, Sex, Reference_Sample
- Optional (reserved) columns: Regions_of_Interest
- See the `sample_table_example.tsv` file (or the sample_table created bye the setup-files command) for a description of individual columns

The config file (default: config.yaml) defines all settings for the analysis and inherits from the inbuilt default.  
Required settings that are not defined by default include array definition files specific to the used array platform and genome build:
- egt_cluster_file: the illumina cluster file (.egt) for the array platform, available from Illumina or the provider running the array 
- bpm_manifest_file: the beadpool manifest file (.bpm) for the array platform, available from Illumina or the provider running the array
- csv_manifest_file (optional): the manifest file in csv format, available from Illumina or the provider running the array

Additionally, the config file needs to define the following paths:
- raw_data_folder: path to the input directory under which the raw data (.idat) can be found. Ths folder should contain subfolders that match the Chip_Name column in the sample table (containing the array chip IDs)
- data_path: the output of StemCNV-check will be written to this path
- log_path: the log files of StemCNV-check will be written to this path


## Usage

Before the first analysis sample table and config file need to be set up (see above). 
Unless otherwise specified, stemcnv-check defaults to look for a "sample_table.tsv" (or .xlsx) and "config.yaml" file.  

Automatic generation of the additional array & genome-build specific static files can only be done if sample data for 
that array is available.  
*Notes*:  

- unless specified directly in the config this will also include download of fasta and gtf file for the reference genome build.
- Array specific files and an updated array_definition block for the config will be written into the cache directory 
  (default: '~/.cache/stemcnv-check'). However, you still need to update or remove the array_definition from your 
  config.yaml file, otherwise the cached definitions and files will not be used.

`stemcnv-check make-staticdata [-s <sample_table>] [-c <config_file>]`

To start the analysis, simply invoke the (implied) run command:

`stemcnv-check [run] [-s <sample_table>] [-c <config_file>]`


## Example data

This repository contains example data (using data from the Genome in a Bottle samples) that can be used to test the setup.
After pulling the repository and creating and activating the base StemCNV-check conda environment, test data can be downloaded via [git LFS](https://git-lfs.com/) and StemCNV-check can be run with the following commands.  
(Note that this will also include the download a fasta and gtf file for the human genome. If you have suitable files available locally, 
it is recommended to replace the corresponding paths in the config.yaml to avoid unnecessary and time-consuming downloads):

Install git lfs and pull test data:
- `sudo apt-get install git-lfs`
- `git lfs fetch`
- `git lfs checkout`

Run the example data:
- `cd example_data`
- `stemcnv-check make-staticdata` 
- `stemcnv-check`

## Output

StemCNV-check will produce the following output files for each sample, when run with default settings:
- `data_path/{sample}/{sample}.annotated-SNP-data.{filter}-filter.vcf.gz`  
  The filtered, processed and annotated SNP data of the array in vcf format
- `data_path/{sample}/{sample}.CNV_calls.CBS.vcf.gz`  
  The CNV calls for the sample from the CBS (Circular Binary Segmentation) algorithm in vcf format
- `data_path/{sample}/{sample}.CNV_calls.PennCNV.vcf.gz`  
  The CNV calls for the sample from the PennCNV caller, in vcf format
- `data_path/{sample}/{sample}.CNV_calls.combined-annotated.vcf.gz`  
  The CNV calls processed, combined and annotated by StemCNV-check, in vcf format. 
  Annotation includes comparison against reference sample; gene annotation; hotspots for stem cells, cancer and dosage 
  sensitivity; call scoring; and call labelling (i.e. as Critical de-novo call).
- `data_path/{sample}/extra_files`  
  Folder in which additional QC log files are stored.
- `data_path/{sample}/{sample}.summary-stats.xlsx`  
  An Excel file with summary information for the sample. The first sheets contains quality summary statistics, including 
  array quality measures, number of CNV and LOH calls, the number of calls above Check_Score thresholds. The further 
  sheets have more details from individual CNV callers or sample comparisons.
- `data_path/{sample}/{sample}.StemCNV-check-report.html`; ... 
  Html report containing summary statistics, QC statistics, lists of CNV calls sorted by annotation score, 
  plots of most/all CNVs and sample comparison. The default 'StemCNV-check-report' only contains plots for the top20 
  calls or calls above a user defined Check_Score threshold. A full report can easily be enabled in the config.yaml. 
  The content of either the default or any additional reports can also be fine-tuned through the config.yaml file.

Furthermore, the following collated summary tables can be created. 
(Optionally with a date prefix, or as tsv instead of xlsx files):
- `data_path/[YYYY-MM-DD_]summary-overview.{xlsx,tsv}`  
  A table that contains the information of the sample wise "summary-stats", but combined for all samples.
  Additionally, more information derived from the sampletable columns can be included. This output is 
  included in the default 'complete' target.
- `data_path/[YYYY-MM-DD_]combined-cnv-calls.{xlsx,tsv}`  
  A table that contains all CNV calls from all samples, that meet config defined filter criteria (By default: 
  everything except calls with a minimum probe/size/density flag). This output is *not* included in the default 
  'complete' target, but can be created with the 'collate-cnv-calls' target.
