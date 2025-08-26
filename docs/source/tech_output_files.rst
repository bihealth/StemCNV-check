.. _basics-output-files:

Output files
============

StemCNV-check will produce the following output files for each sample, when run with default settings:

- ``data_path/{sample}/{sample}.annotated-SNP-data.{filter}-filter.vcf.gz``

  The filtered, processed and annotated SNP data of the array in vcf format

- ``data_path/{sample}/{sample}.CNV_calls.CBS.vcf.gz``

  The CNV calls for the sample from the CBS (Circular Binary Segmentation) algorithm in vcf format

- ``data_path/{sample}/{sample}.CNV_calls.PennCNV.vcf.gz``

  The CNV calls for the sample from the PennCNV caller, in vcf format

- ``data_path/{sample}/{sample}.CNV_calls.combined-annotated.vcf.gz``

  The CNV calls processed, combined and annotated by StemCNV-check, in vcf format. 
  Annotation includes comparison against reference sample; gene annotation; hotspots for stem cells, cancer and dosage 
  sensitivity; call scoring; and call labelling (i.e. as Critical de-novo call).

- ``data_path/{sample}/extra_files``
  
  Folder in which additional QC log files are stored.

- ``data_path/{sample}/{sample}.summary-stats.xlsx``

  An Excel file with summary information for the sample. The first sheets contains quality summary statistics, including 
  array quality measures, number of CNV and LOH calls, the number of calls above Check_Score thresholds. The further 
  sheets have more details from individual CNV callers or sample comparisons.

- ``data_path/{sample}/{sample}.SNV-analysis.xlsx``

  An Excel file with the results from analysis on annotated SNVs (from the SNP probes). This includes a list of all SNVs 
  with an annotated impact on a gene (including the gene name and categorisation based on predicted impact, known hPSC 
  reference SNV hotspots, match and call reliability), coverage of known hotspots by the utilised array, a distance 
  matrix of the sample to other selected samples (based on config settings), and a chromosome based summary of where SNPs occur.

- ``data_path/{sample}/{sample}.StemCNV-check-report.html``; ...

  Html report containing summary statistics, QC statistics, lists of CNV calls sorted by annotation score, 
  plots of most/all CNVs and sample comparison. The default 'StemCNV-check-report' only contains plots for the top20 
  calls or calls above a user defined Check_Score threshold. A fully self-contained report can easily be enabled in the config.yaml. 
  The content of either the default or any additional reports can also be fine-tuned through the config.yaml file.

- ``data_path/{sample}/{sample}.StemCNV-check-report-html_images``

  Folder containing all images included in the html report. 

Furthermore, the following collated summary tables can be created. 
(Optionally with a date prefix, or as tsv instead of xlsx files):

- ``data_path/[YYYY-MM-DD_]summary-overview.{xlsx,tsv}``

  A table that contains the information of the sample wise "summary-stats", but combined for all samples.
  Additionally, more information derived from the sampletable columns can be included. This output is 
  included in the default 'complete' target.

- ``data_path/[YYYY-MM-DD_]combined-cnv-calls.{xlsx,tsv}``

  A table that contains all CNV calls from all samples, that meet config defined filter criteria (By default: 
  everything except calls with a minimum probe/size/density flag). This output is *not* included in the default 
  'complete' target, but can be created with the 'collate-cnv-calls' target.
