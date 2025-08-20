
StemCNV-check Manual
^^^^^^^^^^^^^^^^^^^^

This is the manual for StemCNV-check, a tool written to simplify copy number variation (CNV) analysis of SNP array data, 
specifically for quality control of (pluripotent) stem cell lines.  
StemCNV-check uses snakemake to run the complete analysis from raw data (.idat) up to html report generation for all 
defined samples with a single command. Samples need to be defined in a (tabular) sample table and the workflow settings 
are defined through a yaml config file.


**StemCNV-check features include:**

- Quality control of hPSC genomic integrity based on CNV detection in SNP-array data 
    - Detection of CNVs and loss of heterozygosity
    - Comparison to reference sample
    - Extensive CNV annotation, including ranking by Check-Score 
- Comparison of SNVs against reference and annotation of changes in amino acid sequence
- Comparison of sample identity based on SNP genotypes
    - This allows sample identification and detection of swaps  or cross-contamination
- Easily readable report in html format
    - All results summarised in a single report
    - Summary overview of quality metrics
    - Guided interpretation of results and links to relevant resources

**Requirements:**

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
    :name: stemcnv_tutorials
    
    tut_example_data
    tut_report_analysis
    tut_project_setup    

.. toctree::
    :maxdepth: 2
    :caption: Troubleshooting
    :name: stemcnv_issues

    issues_intro
    issues_faq

.. toctree::
    :maxdepth: 2
    :caption: Technical documentation
    :name: stemcnv_technical
    
    tech_cmd_line
    tech_config
    tech_labels
