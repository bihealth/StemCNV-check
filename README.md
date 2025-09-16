[![CI](https://github.com/bihealth/StemCNV-check/actions/workflows/ci.yml/badge.svg)](https://github.com/bihealth/StemCNV-check/actions/workflows/ci.yml)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/stemcnv-check/README.html)
[![read-the-dcos](https://app.readthedocs.org/projects/stemcnv-check/badge/?version=latest)](https://stemcnv-check.readthedocs.io/)
[![MIT license](https://img.shields.io/badge/License-MIT-green.svg)](https://opensource.org/licenses/MIT)

# StemCNV-check

## About

StemCNV-check is a tool written to simplify copy number variation (CNV) analysis of SNP array data, 
specifically for quality control of (pluripotent) stem cell lines.  

StemCNV-check uses snakemake to run the complete analysis from raw data (.idat) up report generation 
for all defined samples with a single command. Samples need to be defined in a (tabular) sample table and 
the workflow settings are defined through a yaml file.

## Features

StemCNV-check allows generation of read-to-analyse reports from raw SNP-array data in a single command, where possible 
hPSC samples are always compared to a (parental) reference sample and results are annotated with additional information 
to ease interpretation of the data. The report contains and summarises all intermediate analysis steps:

- Summary of quality measures for the sample, including comparison to predefined thresholds
- Sorted list of CNV calls, split by de-novo and reference matching calls
  - Sorting uses the Check-Score from our upcoming manuscript and combines contributions from CNV size and copynumber 
    as well as additions from annotation from overlapping [stem cell hotspots](https://bihealth.github.io/StemCNV-check/CNV-hotspots/index_1.html),
    cancer driver genes, predicted dosage sensitive genes and other gene annotations.
  - The top CNVs have detailed images and tables listing overlapping features and genes with link-outs to further resources
  - Sample specific regions of interest are also displayed in a separate section reargdless of CNV status with full details
- Sorted list of SNVs that affect protein sequences, split by de-novo and reference matching variants
  - Sorting takes known [stem cell hotspots](https://bihealth.github.io/StemCNV-check/SNV-hotspots/index_1.html), sample
    specific regions of interest, SNP quality and severity of predicted protein changes
- Gnome wide overview plots for visual inspection of large aberrations, including side-by-side comparison to reference sample
- Sample identity comparison based on clustering by SNP-genotypes

## Documentation

Please consult our [documentation](https://stemcnv-check.readthedocs.io/) on read-the-docs for detailed instructions on
installation, usage, interpretation, troubleshooting and technical implementation of StemCNV-check.

## Installation

StemCNV-check requires a linux environment (or WSL on windows) and a working conda (or mamba) installation. 
Follow the recommended instructions to install [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) or [mamba](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html).

It is recommended to install StemCNV-check through the bioconda channel. If you do not use conda for other things 
omitting the environment name and installing into your base environment may be an option.

`conda install bioconda::stemcnv-check [-n stemcnv-check]`

## Usage

For detailed explanations on installation, setup and usage please consult our [documentation](https://stemcnv-check.readthedocs.io/).
In brief: 

- Setup your config file and sample table, i.e.: `stemcnv-check setup-files`
- Once only, use the `stemcnv-check make-staticdata` to create static data for your specific array
- Start the analysis with the run command: `stemcnv-check run`

### Example data (legacy)

We provide example data, so you can easily test StemCNV-check for yourself, please consult the documentation for details. 
This data is primarily available through [zenodo](https://zenodo.org/records/16962381), 
however a legacy version through git-lfs also still exists.

Install git lfs and pull test data:
- `sudo apt-get install git-lfs`
- `git lfs fetch`
- `git lfs checkout`

Run the example data:
- `cd example_data`
- `stemcnv-check make-staticdata` 
- `stemcnv-check run`



