

 QC measure by stemCNV for the array
colour codes in the report (3 categories: green, yellow, red)
measures:
Call rate
SNP_pairwise distance to reference: [500, 5000]
-measured as distance btw reference and sample
-uses minimum threshold
-depends on array platform

Loss_gain_log2ratio: [2, 4]
Total calls CNVs
Total calls LOH
Critical calls CNV
Critical cals LOH
Critical CNVs

**Specification of labels (and their report colors) assigned to sample level QC measures**
sample_labels:
    OK: green
    unusual: yellow
    warning: orange
    high concern: red

**Possible/Defined report sections**
report_sections:
  - sample.information
  - QC.summary
  - QC.GenCall
  - QC.PennCNV
  - QC.CBS
  - QC.settings
  - SNV.table
  - SNV.hotspot.coverage
  - SNV.QC.details
  - denovo_calls.table
  - denovo_calls.plots
  - reference_gt_calls.table
  - reference_gt_calls.plots
  - regions.of.interest
  - SNP.dendrogram
  - genome.overview

**Possible/Defined subsections in the CNV plot sections**
report_plotsections:
  - denovo
  - reference_gt
  - regions_of_interest


**Possible/Defined QC measures on sample level**
sample_qc_measures:
    - call_rate
    - computed_gender
    - SNPs_post_filter
    - SNP_pairwise_distance_to_reference
    - loss_gain_log2ratio
    - total_calls_CNV
    - total_calls_LOH
    - reportable_calls_CNV
    - reportable_calls_LOH
    - reportable_SNVs
    - critical_calls_CNV
    - critical_calls_LOH
    - critical_SNVs