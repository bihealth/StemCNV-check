StemCNV-check overview 
=============

StemCNV-check is a tool written to simplify copy number variation (CNV) analysis of SNP array data, specifically for quality control of (pluripotent) stem cell lines. StemCNV-check uses snakemake to run the complete analysis from raw data (.idat) up report generation for all defined samples with a single command. Samples need to be defined in a (tabular) sample table and the workflow settings are defined through a yaml file.


Workflow  
=============

- Raw data processing 
    * idat files
      
- CNV calling and processing 
- CNV scoring and evaluation 
- SNP & SNV analysis 
- Report generation 
    
.. image:: workflow.png
   :width: 600
