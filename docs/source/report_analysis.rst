Report analysis 
===========================

Report sections       
----------------

- Sample Overview 


   - Sample information
   - QC metrics 
- CNV calling
- SNV analysis
- Sample comparison

1. Sample overview 
----------------

.. image:: sample_info.png
   :width: 600


**Sample information table**

- Contains information from the sample table about the sample named by sample_id: Sex, Reference_Sample, Array_Name, Chip_Name, Chip_Pos

- Reference sample is the sample for comparison, it is a precursor cell line or earliest progenitor cell line for this sample with iPSC clone. It is defined by the sample id in the sample table. 

.. image:: sample_table.png
   :width: 600
The sample table with all samples specifies the sample id and the reference for analysis of each sample. A sample table is used as input for running StemCNV-check.

**QC metrics**

- Summary (Data and Sample QC explanations)

- GenCall
- PennCNV
- CBS
- StemCNV-check utilizes both the PennCNV and CBS algorithms
- Config
- R session info

**Summary.** Data and Sample QC explanations: these summary tables are meant to serve as a quick overview of the quality of an hPSC sample. 

**Data QC explanations:** QC metrics primarily related to the SNP data quality (affected by both the DNA used and the array run itself), this table will also display values from the reference sample if possible. 

**Sample QC explanations:** QC metrics related to the potentially problematic CNVs and SNVs identified in only the analysed sample. This table sums up all variant findings from the analysed sample, which were flagged as critical or reportable.

Note that in contrast to general SNP probes on the array, only those single variants that actually show an alternative allele and affect a protein are considered SNVs by StemCNV-check. Variants that match the genotype of assigned reference samples are never considered critical or reportable.




Data QC measures table
===========================

.. image:: qc_metrics.png
   :width: 600

.. image:: coloring.png
   :width: 500

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

.. image:: data_qc.png
   :width: 600

.. image:: sample_qc.png
   :width: 600

- **Call rate** is % of loci (SNP, CNV) genotyped for the sample. Call rate > 0.99 (default threshold), indicates good-quality data.

.. math:: Call rate = \frac{called markers}{all markers}

For high-quality data 99.5% call rate is expected. However, accuracy is highly sample dependent. When samples do not perform as expected, experimenters can choose to reprocess these samples to confirm or potentially improve results. Poorly performing samples can be systematically excluded from the project. 

- **Computed gender:** M (male) or F (female), should match the value in “Sex” column from the sample table;

- **SNPs Post Filter:** “good quality” SNPs that passed the QC thresholds;

- **SNP Pairwise distance to reference:** absolute GT distance between a sample and its reference. It reflects the similarity between the two cell lines. The smaller the distance (number of different SNPs) the smaller the phylogenetic distance (higher genetic relation between the samples).
- **Loss Gain Log2ratio:** difference in SNP signal intensity between the sample and the reference
- **Total calls CNV:** number of CNVs detected 
- **Total calls LOH:** number of LOH regions detected 

CNVs (copy number variation) are increases or decreases in chromosomal copies for a given region in the genome;

LOH (loss of heterozygosity): a region that no longer has two different alleles has a LOH;
Homozygosity: a locus can duplicate one chromosome and transpose it to the other chromosome;
Hemizygosity: a region can be deleted entirely, leaving only one chromosomal copy;




