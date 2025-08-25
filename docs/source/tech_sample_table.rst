.. _tech-sample-table:

Sample table options
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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
------------------------------------------

- **Required Columns**: Sample_ID, Chip_Name, Chip_Pos, Array_Name, Sex, Reference_Sample, Regions_of_Interest, Sample_Group

Specific explanations for columns:

  - Sample_ID

    The folder and file names for samples are derived from this entry. All entries *must* be unique. 
    To prevent issues with filenames only alphanumeric characters (all letters and number) and the characters ``-`` 
    and ``_``(dash and underscore) are allowed by default.

  - Chip_Name and Chip_Pos

    These entries must match the Sentrix name (usually a 12 digit number) and position (usually ``R..C..``) on the Illumina array

  - Array_Name

    The name of the array used for the sample. This needs to match one of the arrays defined in the config under ``array_definition``

  - Sex

    The sex of the sample is needed for analysis and mandatory. Allowed values are: ``f``, ``female``, ``m`` and ``male`` (not case sensitive)

  - Reference_Sample

    This column should refer to the (exact) Sample_ID of reference sample (i.e. a parental fibroblast line or master bank),
    if there is no usable or applicable reference sample the entry should be empty (i.e. for fibroblast samples).   
    Reference samples are assumed to be the clonal "parents" of a sample. 

  - Regions_of_Interest

    This column can define sample specific regions of interest (i.e. gene edited sites), if none exist ic can be left 
    empty but still must be included in the samplesheet.  
    The syntax for regions of interest is ``NAME|region``, the ``NAME|`` part is optional and mainly useful for 
    labeling or describing the region.  
    The ``region`` part is mandatory and can be one of the following:  
    1) Position: "chrN:start-end": ``chrN`` can be i.e. ``chr3`` or just ``3``,
       start and end are coordinates (which are genome build specific!)
    2) Genomic band: i.e. "4q21.3": a cytogenetic band, both full bands (q21) and subbands (q21.3) are allowed 
    3) Gene symbol: i.e. "TP53": The gene name (or symbol) needs to exactly match the reference annotation (UCSC gtf)
    
    Multiple regions for a single sample should all be in one column entry and be separated by a ``;``

  - Sample_Group

    This column can be used for annotation of similar samples.  
    By default all samples with the same entry will be included in sample comparison based on SNP clustering.

								
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
