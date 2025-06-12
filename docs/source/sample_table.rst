Setting up the sample table
============

The sample table (default: sample_table.tsv) is a tab-separated file describing all samples to be analyzed.
**Excel or tsv** formats are supported.

Empty example files for the sample table and config can be created with this command:

``stemcnv-check setup-files``

If you prefer to use an xlsx file here you can create an example by using:

``stemcnv-check setup-files --sampletable-format xlsx``

You can also use your own Excel file, if the following criteria are met:

  - The actual sample table is in the first sheet of the file and this sheet *only* contains columns for the sample table (optionally with commented lines starting with a '#')

  - All required columns are present and correctly named (the order of columns is not important)
  - It is possible to deviate from the standard column names, but the expected column names need be contained in the actual column names and there needs to a singular way to extract them (via regex).
  - In this case you need to use the ``--column-remove-regex`` option to tell the pipeline how to modify your column names to derive the expected names. If used without an explicit regex (for expert users) spaces and anything following them will be removed from your column names.

  - A simple example with ``--column-remove-regex`` (default) option would be to use i.e:
    'Sample_ID for pipeline', 'Chip_Name (Sentrix Barcode)', 'Chip_Pos (Sentrix Position)'

Filling in the sample table with your data
----------

- **Required Columns**: Sample_ID, Chip_Name, Chip_Pos, Array_Name, Sex, Reference_Sample, Regions_of_Interest, Sample_Group

Specific explanations for columns:
 - Sample_ID:
       The folder and samples names for samples are derived from this entry. All entries *must* be unique.
       To prevent issues with filenames only alphanumeric characters (all letters and number) and the characters -_
       (dash and underscore) are allowed. Include bank ID when possible, only: - or _, do not use special characters: (), {}, /, \, ~,*, & Name has to be UNIQUE.
       This column has auto-formatting enabled, so that the IDs will work with the CNV-pipeline:

       - red entries are either duplicate or contain not-allowed characters (/ and .\)

       - orange entries contain characters that the pipeline will remove (since they can cause issues if used in file names):  :,;()[]{}!?* and <space>
 - Chip_Name and Chip_Pos:
       These entries must match the Sentrix name (usually a 12 digit number) and position (usually R..C..) on the Illumina array
 - Array_Name
       The name of the array used for the sample. This needs to match one of the arrays defined in the config under `array_definition`
 - Sex
       The sex of the sample is needed for analysis and mandatory. Allowed: f[emale]/m[ale] (not case sensitive)
 - Reference_Sample
       This column should refer to the (exact) Sample_ID of reference sample (i.e. a parental fibroblast line or master bank)
      If there is no usable or applicable reference sample the entry should be empty
 - Regions_of_Interest
       Definition of regions for which plots are always generated in the report (i.e. gene edited sites)
       The syntax for regions of interest is `NAME|region`, the `NAME|` part is optional and mainly useful for
       labeling or describing the region.
       The `{region}` part is mandatory and can be one of the following:
       1) Position, "chrN:start-end": `chrN` can be i.e. 'chr3' or just '3', start and end are coordinates (which are genome build specific!)
       2) Genomic band, i.e. "4q21.3": a cytogenetic band, both full bands (q21) and subbands (q21.3) are allowed
       3) Gene symbol, i.e. "TP53": The gene name (or symbol) needs to exactly match the reference annotation (UCSC gtf)
       Multiple regions for a single sample should all be in one column entry and be separated by a `;`
 - Sample_Group
       This column can be used for annotation samples is used by default to select samples for clustering by SNPs.


								
.. list-table::  Example Sample table
   :widths: 15 15 10 10 10 10 10 10 10 
   :header-rows: 1
								
   * - Sample_ID 
     - Chip_Name
     - Chip_Pos
     - Array_Name
     - Sex
     - Reference_Sample
     - Regions_of_Interest
     - Sample_Group
     - Coriell_ID
   * - HG001
     - 207521920117
     - R09C02
     - ExampleArray
     - female
     -
     -
     - 
     - NA12878
   * - HG002
     - 207521920117
     - R05C02
     - ExampleArray
     - male
     -
     -
     - 
     - NA24385
   * - HG004
     - 207521920117
     - R07C02
     - ExampleArray
     - female				
     -
     -
     - 
     - NA24143
   * - HG005
     - 207521920117
     - R01C02
     - ExampleArray
     - male
     -
     -
     - HG006
     - NA24631
   * - HG006
     - 207521920117
     - R03C02
     - ExampleArray
     - male
     -
     -
     - 
     - NA24694
   * - HG007
     - 207521920117
     - R11C02
     - ExampleArray
     - female
     -
     -
     - 
     - NA24695

**Extended sample table. Additional data types/columns.**

- Any number of additional columns can be added to the sample table as well, unless referred to in the config they will be ignored.

- Line family (iPSC line names without the clone part)	
- DNA ID/ Barcode (CORE)	
- Gender	
- Passage	
- Gene edited (yes/no)	
- Passages after editing	
- Type of editing	
- `Modification <https://scc-docs.charite.de/openkm/kcenter/#/browser/uuid/6f505d68-4e61-4f2d-a46d-4ad434ea94d5>`_ . Check Gene Editing Overview table to input correct modification
- Chromosome	
- ROI for StemCNV-Check	
- Bank	(Only use: MBXX WBXX seed primary)
- Cell type (iPSC/reference)
- latest parental CONTROL sample (patient cells or preceeding Bank MB/WB/Seed). If it is not 'reference' then sample name chosen for this column MUST exist in the first column
- earliest parental CONTROL (patient cells or MB). If it is not 'reference' then sample name chosen for this column MUST exist in the first column
- AG (resp user)	
- Service request ID openIRIS	
- Responsible person (CORE)	
- Batch group	
- Additional references (e.g. for dendrogram). This column works the same as the "Parental Control" one, except that you can add multiple references separated by commas (in the same field). Excel can not do conditional formatting for that.
- Send to L&B (date)	
- Data received (date)	
- Sample_Name (L&B)	
- Chip/Sentrix Barcode (L&B)	
- SentrixPosition (L&B)	
- Chip Type (L&B)	
- Manifest Version	
- Pass/fail (Use pass/fail ONLY for non-reference samples!!)
- Analysis by	
- Report generated/  updated	
- Results/Comment	
- known CNVs in this line	
- Sample derived from	
- Culture medium (used for routine maintenance culture)	Coating	Hypoxya (5% O2)/ Normoxya (20% O2)	
- Passaging method (for routine maintenance)	
- Survival factor for enzymatic passaging (maintenance)	

- Reprogramming method
