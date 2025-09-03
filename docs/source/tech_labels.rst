Assigned labels and filters
===========================

.. caution::Under construction
    This page is still under construction and has not been finalised yet

.. _tech-cnv-labels:

CNV labels
----------

CNV calls in the vcf files are annotated with potential *FILTER* values, 
that describe a specific reason a call may not be trustworthy.

.. list-table::
   :widths: 20 30
   :header-rows: 1

   * - Filter name
     - Description

   * - min_size
     - CNV call below minimum size (<1000bp)
   * - min_probes
     - CNV call from <5 probes
   * - min_density 
     - CNV call with <10 probes/Mb
   * - high_probe_dens
     - Probe density of segment is higher than 99% of the array
   * - probe_gap
     - Probe coverage of segment has considerable gap (min. 33% depending on probe number - see config)

All CNV calls are given a label based on their check score, filters and reference match.

- The labels described here are always available, but can be changed or new labels can be added through the config file
- If not other category fits (which should not occur with default settings), then the last defined **"Exclude call"** label will always be assigned.

.. list-table::  
   :widths: 18 10 25 10 25
   :header-rows: 1

   * - CNV label
     - Minimum Check_Score
     - FILTER values not allowed
     - Match in Reference Sample
     - Description

   * - Critical de-novo
     - ≥55
     - high_probe_dens, probe_gap, min_size, min_probes, min_density
     - No
     - High-confidence CNV call with likely biological relevance
   * - Reportable de-novo
     - ≥55
     - min_size, min_probes, min_density
     - No
     - CNV call with potential biological relevance
   * - de-novo call
     - ≥ 0
     - min_size, min_probes, min_density
     - No
     - General CNV call meeting minimal quality requirements
   * - Reference genotype
     - ≥ 0 (any)
     -
     - Yes
     - CNV call that matches the reference sample genotype
   * - Excluded call
     - ≥ 0 (any)
     -
     - No
     - CNV call that does not match minimal quality requirements (close to "noise")

.. _tech-snv-labels:

SNV labels
----------

Each SNV is first assigned to a category based on annotation and overlap with hPSC hotspots or sample specific ROIs>

.. list-table::
   :widths: 20 60
   :header-rows: 1

   * - SNV_category
     - Description

   * - ROI-overlap
     - SNV overlapping a sample specific regions of interest

   * - hotspot-match
     - SNV matching a known stem cell hotspot mutation (see also SNV hotspot coverage)

   * - hotspot-gene
     - SNV overlapping a sample specific regions of interest

   * - protein-ablation
     - SNV (likely) fully disrupting protein function (i.e. frameshift, stop gain, stop loss)

   * - protein-changing
     - SNV causing a change the protein sequence (i.e. missense, inframe)

   * - other
     - SNV with other unclear or undetermined effect on protein function


SNVs are then assigned a label based on their category, reference match and genotype quality:

.. list-table::  SNV labels
   :widths: 12 12 12 12 50
   :header-rows: 1

   * - SNV Label
     - Categories     
     - Match in Reference Sample
     - Genotype quality
   * - Critical
     - hotspot-match
     - No
     - high
   * - Reportable
     - ROI-overlap, hotspot-gene, protein-ablation
     - No
     - high
   * - Unreliable critical/reportable
     - hotspot-match, ROI-overlap, hotspot-gene, protein-ablation
     - No
     - low
   * - de-Novo
     - protein-changing, other
     - No
     - any
   * - Reference genotype
     - any
     - Yes
     - any

