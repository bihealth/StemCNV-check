
StemCNV-check Manual
^^^^^^^^^^^^^^^^^^^^

.. admonition:: Under construction

    This documentation is still under construction and not all parts have been finalised yet

This is the manual for StemCNV-check, a tool written to simplify copy number variation (CNV) analysis of SNP array data, 
specifically for quality control of (pluripotent) stem cell lines. 

StemCNV-check uses snakemake to run the complete analysis from raw data (.idat) up to html report generation for all 
defined samples with a single command. Before this one command, samples need to be defined in a (tabular) sample table, 
the workflow settings through a yaml config file and array specific static data needs to be generated once.


**StemCNV-check features include:**

- Quality control of hPSC genomic integrity based on CNV detection in SNP-array data

  - Detection of CNVs and loss of heterozygosity
  - Comparison to reference sample
  - Extensive CNV annotation, including ranking by Check-Score

- Comparison of SNVs against reference and annotation of changes in amino acid sequence
- Comparison of sample identity based on SNP genotypes
- Easily readable report in html format

  - All results summarised in a single report
  - Summary overview of quality metrics
  - Guided interpretation of results and links to relevant resources


Requirements
------------

StemCNV-check requires a linux environment due to dependencies on open source packages, 
for windows users WSL (Windows Subsystem for Linux) is recommended. 
Easy installation is possible through the bioconda repository. 
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
    basics_output_files

 
.. toctree::
    :maxdepth: 1
    :caption: Tutorials
    :name: stemcnv_tutorial

    tut_example_data
    tut_report_analysis
    tut_project_setup    

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
    tech_cmd_line
    tech_cache
    tech_config
    tech_labels
