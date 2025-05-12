
**Evaluation settings**
---------------

All CNV calls are given a label based on their check score, filters and reference match. The labels described here are always available, but can be changed or new labels can be added. If not other category fits (which should not occur with default settings), then the last defined "Exclude call" label will always be assigned.

  •  Possible values for the "not_allowed_vcf_filters" list are: **probe_gap, high_probe_dens, min_size, min_probes, min_density**

.. list-table::  CNV_call_labels
   :widths: 25 25 25 25
   :header-rows: 1

   * - CNV_call_labels
     - minimum_check_score
     - not_allowed_vcf_filters
     - reference_match

   * - Critical de-novo
     - 55
     -  ['high_probe_dens', 'probe_gap', 'min_size', 'min_probes', 'min_density']
     - FALSE
   * - Reportable de-novo
     - 55
     -  ['min_size', 'min_probes', 'min_density']
     - FALSE
   * - de-novo call
     - 0
     - ['min_size', 'min_probes', 'min_density']
     - FALSE
   * - Reference genotype
     - 0
     - []
     - TRUE
   * - Excluded call
     - 0
     - []
     - FALSE


**Labelling system**

**Specification of labels (and their report colors) assigned to sample level QC measures**
sample_labels:
    OK: green
    unusual: yellow
    warning: orange
    high concern: red

# Default labels for CNVs (more can be added by users)
CNV_labels:
    # This is used to count critical CNVs & LOHs
    - Critical de-novo
    # This is used to count reportable CNVs & LOHs
    - Reportable de-novo
    - de-novo call
    - Reference genotype
    - Excluded call

# possible labels for SNVs
SNV_labels:
    - critical
    - reportable
    - unreliable impact
    - de-novo SNV
    - reference genotype

**Label for CNVs merged from multiple callers**

combined_cnvs: 'combined-call'

**The following lists are primarily used by the check_config functions**

Possible/Defined FILTERs applied to CNV calls (vcf style)

vcf_filters:
    - probe_gap
    - high_probe_dens
    - min_size
    - min_probes
    - min_density

**Possible/Defined categories for SNVs, each category can be assigned critical or reportable**
SNV_category_labels:
    - ROI-overlap
    - hotspot-match
    - hotspot-gene
    - protein-ablation
    - protein-changing
    - other

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