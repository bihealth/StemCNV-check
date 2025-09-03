Assigned labels and filters
===========================

.. caution::Under construction
    This page is still under construction and has not been finalised yet

.. _tech-cnv-labels:

CNV labels
----------

**CNV_labels:**
    - Critical de-novo (used to count critical CNVs & LOHs)
    - Reportable de-novo (used to count reportable CNVs & LOHs)
    - de-novo call
    - Reference genotype
    - Excluded call

Possible/Defined filters applied to CNV calls (vcf style) **vcf_filters:**
  - probe_gap
  - high_probe_dens
  - min_size
  - min_probes
  - min_density

- All CNV calls are given a label based on their check score, filters and reference match.
- The labels described here are always available, but can be changed or new labels can be added.
- If not other category fits (which should not occur with default settings), then the last defined **"Exclude call"** label will always be assigned.
- Possible values for the Exclusion filters applied ("not_allowed_vcf_filters"): **probe_gap, high_probe_dens, min_size, min_probes, min_density**

.. list-table::  CNV_call_labels
   :widths: 18 10 25 10 10 25
   :header-rows: 1

   * - CNV_call_labels
     - Minimum Check_Score
     - Exclusion Filters applied
     - Match in Reference Sample
     - Interpretation/ Concern
     - Description

   * - Critical de-novo
     - ≥55
     - high_probe_dens, probe_gap, min_size, min_probes, min_density
     - No
     - High
     - High-confidence CNV call with biological relevance
   * - Reportable de-novo
     - ≥55
     - min_size, min_probes, min_density
     - No
     - Middle
     - CNV call with potential biological relevance
   * - de-novo call
     - ≥ 0
     - min_size, min_probes, min_density
     - No
     - Low
     - General CNV call meeting minimal quality requirements
   * - Reference genotype
     - ≥ 0 (any)
     -
     - Yes
     - None
     - CNV call that matches the reference sample genotype
   * - Excluded call
     - ≥ 0 (any)
     -
     - No
     - None
     - CNV call that does not match minimal quality requirements (close to "noise")

**Label for CNVs merged from multiple callers** - combined_cnvs: 'combined-call'

.. list-table::  Exclusion filter
   :widths: 20 30 10 10 30
   :header-rows: 1

   * - Exclusion filter
     - Description
     - Basic
     - Extended
     - Description

   * - min_size (minimal size of CNV)
     - <1000bp
     - X
     - X
     - CNV call below minimum size (<1000bp)
   * - min_probes (minimal number of probes )
     - <5 probes
     - X
     - X
     - CNV call from <5 probes
   * - min_density (minimal probe density per call size)
     - <10 per 1Mb
     - X
     - X
     - CNV call with <10 probes/Mb
   * - high_probe_dens (too high probe density per call size)
     - density >99% than the array
     -
     - X
     - Probe density of segment is higher than 99% of the array
   * - probe_gap (presence of a gap between probes)
     - no probe region between two probes (min. 33% of call, more with low probe count)
     -
     - X
     - Probe coverage of segment has considerable gap (min. 33% depending on probe number - see config)

.. _tech-snv-labels:

SNV labels
----------

**SNV_labels**:

- **critical**
                 SNV with likely critical significance on hiPSC line

- **reportable**
                 reportable	SNV with possible significance on hiPSC line
- **unreliable impact**
                 unreliable impact	SNV with likely or possible significance on hiPSC line, but unreliable signal
- **de-novo SNV**
                  SNV with de-novo status, but no clear functional impact
- reference genotype
                 SNV already detected in the reference sample


**Possible/Defined categories for SNVs**

.. list-table::  SNV categories
   :widths: 20 60 20
   :header-rows: 1

   * - SNV_category_labels
     - Description
     - “Critical” / “reportable”

   * - ROI-overlap
     - SNV overlapping a sample specific regions of interest
     -                 Critical
   * - hotspot-match
     - SNV matching a known stem cell hotspot mutation (see also SNV hotspot coverage)
     -                 Critical
   * - hotspot-gene
     - SNV overlapping a sample specific regions of interest
     -                 Reportable
   * - protein-ablation
     - SNV (likely) fully disrupting protein function (i.e. frameshift, stop gain, stop loss)
     -                 Reportable
   * - protein-changing
     - SNV causing a change the protein sequence (i.e. missense, inframe)
     -
   * - other
     - SNV with other unclear or undetermined effect on protein function
     -


Each category can be assigned critical or reportable.

.. list-table::  SNV labels
   :widths: 12 12 12 12 50
   :header-rows: 1

   * - SNV Label
     - Match in Reference Sample
     - Impact
     - Interpretation/ Concern
     - Description

   * - Critical
     - No
     - High/moderate
     - High/ middle
     - Overlaps a specific ROI or matches/overlap a described stem cell hotspot mutation. Could affect protein function.
   * - Reportable
     - No
     - High/moderate
     - Middle
     - **High**- causes protein loss of function in any gene (incl. gene with SNV-stem cell hotspot). **Moderate**-causes protein change in a gene that also has a SNV-stem cell hotspot.
   * - Unreliable critical/reportable
     - No
     - High/moderate
     - Low
     - Could belong to critical or reportable categories but technical scores for genotype are low e.g. bad quality calls or could be missing in the reference. If of concern, genotype should be confirmed by another method.
   * - de-Novo
     - No
     - High/moderate
     - Low
     - Causes protein change but biological impact is unknown
   * - Reference genotype
     - Yes
     - High/moderate
     - None
     - Also detected in the reference sample

