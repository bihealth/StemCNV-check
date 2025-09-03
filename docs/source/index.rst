.. _index:

StemCNV-check Manual
^^^^^^^^^^^^^^^^^^^^

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

.. admonition:: Under construction

    This documentation is still under construction and not all parts have been finalised yet


Requirements
------------

StemCNV-check can run on Linux and Windows systems. On Windows it requires the Windows Subsystem for Linux (WSL) due to dependencies on open source packages. Easy installation is possible through the bioconda repository. 
We are working towards native support for all major operating systems for our bioconda release.

StemCNV-check requires at least 8GB of RAM and 4 CPU cores with default settings, however 12-16GB RAM are recommended.
Running StemCNV-check with fewer resources may work, but has not been extensively tested.

Table of Contents
^^^^^^^^^^^^^^^^^

* :ref:`index`

.. toctree::
    :maxdepth: 1
    :caption: Introduction

    intro_workflow
    intro_quickstart

.. toctree::
    :maxdepth: 2
    :caption: Step-by-step instructions
    :name: stemcnv_install

    steps1_nonlinux
    steps2_installation
    steps3_basic_usage
    steps4_example_data
    steps5_report_analysis
    steps_ext_project_setup    

.. toctree::
    :maxdepth: 1
    :caption: Troubleshooting
    :name: stemcnv_issues

    issues_intro
    issues_common
    issues_support

.. toctree::
    :maxdepth: 1
    :caption: Technical documentation
    :name: stemcnv_technical

    tech_output_files
    tech_cache    
    tech_cmd_line
    tech_config
    tech_labels
