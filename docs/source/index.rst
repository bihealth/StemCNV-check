
StemCNV-check Manual
^^^^^^^^^^^^^^^^^^^^

.. admonition:: Under construction

    This documentation is still under construction and not all parts have been finalised yet

This is the manual for StemCNV-check, a tool developed to simplify copy number variation (CNV) analysis using SNP array data, 
specifically for quality control of human pluripotent stem cells (hPSC). 

**StemCNV-check features include:**

- Quality control of hPSC genomic integrity based on CNV detection in SNP-array data

  - Detection of CNVs and loss of heterozygosity (LOH)
  - Comparison to a reference sample
  - Extensive CNV annotation, including ranking by Check-Score

- Comparison of single nuleotid variants (SNV) against a reference and annotation of changes in amino acid sequence
- Comparison of sample identity based on SNP genotypes
- Easily readable report in html format

  - All results summarised in a single report
  - Summary overview of quality metrics
  - Guided interpretation of results and links to relevant resources

StemCNV-check uses snakemake to run the complete analysis from raw data (.idat) up to html report generation for all 
defined samples with a single command. Before this one command, samples need to be defined in a (tabular) sample table, 
the workflow settings through a yaml config file and array specific static data needs to be generated once.

Requirements
------------

StemCNV-check can run on Linux and Windows systems. On Windows it requires the Windows Subsystem for Linux (WSL) due to dependencies on open source packages. Easy installation is possible through the bioconda repository. 
We are working towards native support for all major operating systems for our bioconda release.


Table of Contents
^^^^^^^^^^^^^^^^^

.. toctree::
    :maxdepth: 2
    :caption: Contents:

.. toctree::
    :maxdepth: 2
    :caption: Getting started 
    :name: stemcnv_install

    basics_nonlinux
    basics_installation
    basics_usage

.. toctree::
    :maxdepth: 1
    :caption: Tutorials
    :name: stemcnv_tutorial

    tut_example_data
    tut_project_setup 
    tut_report_analysis
   
.. toctree::
    :maxdepth: 2
    :caption: Troubleshooting
    :name: stemcnv_issues

    issues_intro
    issues_common
    issues_support

.. toctree::
    :maxdepth: 2
    :caption: Technical documentation
    :name: stemcnv_technical
    
    tech_workflow
    basics_output_files
    tech_cmd_line
    tech_cache
    tech_config
    tech_labels
