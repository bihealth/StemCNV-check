
===========================
Report sections       
===========================

- Sample Overview 


   - Sample information
   - QC metrics 
- CNV calling
- SNV analysis
- Sample comparison

1. Sample overview 
===========================

.. image:: sample_info.png
   :width: 600


**Sample information table**

- Contains information from the sample table about the sample named by sample_id: Sex, Reference_Sample, Array_Name, Chip_Name, Chip_Pos

- Reference sample is the sample for comparison, it is a precursor cell line or earliest progenitor cell line for this sample with iPSC clone. It is defined by the sample id in the sample table. 

.. image:: sample_table.png
   :width: 600
The sample table with all samples specifies the sample id and the reference for analysis of each sample. A sample table is used as input for running StemCNV-check.


QC measures 
===========================

**Sections** 

- **Summary** with two tables (Data and Sample QC explanation)Data and Sample QC explanations: these summary tables are meant to serve as a quick overview of the quality of an hPSC sample. 

**Data QC explanations:** QC metrics primarily related to the SNP data quality (affected by both the DNA used and the array run itself), this table will also display values from the reference sample if possible. 

**Sample QC explanations:** QC metrics related to the potentially problematic CNVs and SNVs identified in only the analysed sample. This table sums up all variant findings from the analysed sample, which were flagged as critical or reportable.

Note that in contrast to general SNP probes on the array, only those single variants that actually show an alternative allele and affect a protein are considered SNVs by StemCNV-check. Variants that match the genotype of assigned reference samples are never considered critical or reportable.

- **GenCall** (Illumina genotyping algorythm)

- **PennCNV and CBS** (two CNV calling algorithms utilized by StemCNV-check)
- **Config** (file with analysis parameters)

- **R session info** (output from terminal in R, useful for determining errors)


.. image:: qc_metrics.png
   :width: 600

.. image:: coloring.png
   :width: 500

Data QC explanation 
-------- 

.. image:: data_qc.png
   :width: 600


- **Call rate** is % of loci (SNP, CNV) genotyped for the sample. Call rate > 0.99 (default threshold), indicates good-quality data.

.. role:: raw-math(raw)
    :format: latex html

:raw-math:`$$ \text{Call rate} = \frac{\text{called markers}}{\text{all markers}} $$`

For high-quality data 99.5% call rate is expected. However, accuracy is highly sample dependent. When samples do not perform as expected, experimenters can choose to reprocess these samples to confirm or potentially improve results. Poorly performing samples can be systematically excluded from the project. 

- **Computed gender:** M (male) or F (female), should match the value in “Sex” column from the sample table;

- **SNPs Post Filter:** “good quality” SNPs that passed the QC thresholds;

- **SNP Pairwise distance to reference:** absolute GT distance between a sample and its reference. It reflects the similarity between the two cell lines. The smaller the distance (number of different SNPs) the smaller the phylogenetic distance (higher genetic relation between the samples).

- **Total calls CNV:** number of CNVs detected.
 | **CNVs** (copy number variation) are increases or decreases in chromosomal copies for a given region in the genome;
- **Total calls LOH:** number of LOH regions detected 
 | **LOH (loss of heterozygosity):** a region that no longer has two different alleles has a LOH;
 | **Homozygosity:** a locus can duplicate one chromosome and transpose it to the other chromosome;
 | **Hemizygosity:** a region can be deleted entirely, leaving only one chromosomal copy;

- **Loss Gain Log2 ratio** (signal intensity - log R ratio): difference in SNP signal intensity between the sample and the reference
 | Log R Ratio is a metric that normalises signal intensity for CNV analysis. It represents the number of copies relative to the normal reference sample.  
 | LRR deviation from an average of 0 **indicates a gain or a loss.**

.. role:: raw-math(raw)
    :format: latex html

:raw-math:`$$ \text{log R Ratio} = \log{2}{\frac{\text{CNV signal  intensity sample}}{\text{CNV signal intensity  ref}}} $$`

LRR deviation from an average of 0 indicates a gain or a loss.

.. list-table::  
   :widths: 50 50
   :header-rows: 0

   * - Log2 Ratio = 0
     - SNP has the expected copy number (usually two in a diploid genome)
   
   * - Log2 Ratio > 0
     - Indicates a gain in copy number, meaning the sample has more than two copies of the region
 
   * - Log2 Ratio < 0
     - Indicates a loss in copy number, meaning the sample has fewer than two copies of the region


**Threshold measures set in the config file (can be changed by user):**

- **call rate**: [0.99, 0.99]

- **SNP_pairwise_distance_to_reference:** [500, 5000]. It is based on the array platform. [500,5000] for GSA array (~700k probes).
- **loss_gain_log2ratio:** [2, 4]
- **total_calls_CNV:** [10, 50]
- **total_calls_LOH:** [30, 75]
- **reportable_calls_CNV:** [5, 10]
- **reportable_calls_LOH:** [5, 10]
- **critical_calls_CNV:** [1, 1]
- **critical_calls_LOH:** [1, 1]
- **reportable_SNVs:** [5, 10]

- **critical_SNVs:** [1, 1]


Sample QC explanation  
-------- 

- Reportable calls CNV

- Reportable calls LOH 
- Reportable SNVs
- Critical calls CNV
- Critical calls LOH 

- Critical SNVs

.. image:: sample_qc.png
   :width: 600

| **Call frequency**: % of samples that a genotype was called for the SNP

BAF and log2 ratio charts
--------
| **B allele frequency (BAF)** is the proportion of the B allele signal relative to the total signal for a SNP. It is  a normalized representation of how often B allele is called.

- **Normal, heterozygous samples**: three distinct bands are seen. Homozygous calls are at the top (1.0) and bottom (0.0) of the chart, representing the BB and AA calls, respectively. The middle band at 0.5 represents AB, 50% of BAF in the genotype. 
- **Loss of heterozygosity (LOH)**: the middle band is missing while the bands at 1.0 and 0.0 remain. BAF of 1.0 can mean either a homozygous genotype of BB or a hemizygous genotype of B [-].
| *BAF charts alone can not distinguish copy neutral LOH from deletion events. That’s why it is necessary to look at the Log ratio chart simultaneously.*
