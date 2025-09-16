.. _steps-ext-project-setup:

Extended project setup
^^^^^^^^^^^^^^^^^^^^^^

Before the first analysis sample table and config file need to be set up. Unless otherwise specified, StemCNV-check 
defaults to look for a ``sample_table.tsv`` (or ``sample_table.xlsx``) and `config.yaml` file. 
To use sample table and config files with different names different (i.e. to compare the results of different settings), 
the command line options ``--sample-table`` or ``--config`` can be used.

It is recommended to start by **creating a separate folder** for your project. This folder should include raw data folder, config.yaml and sample table files.

Empty example files for the sample table and config can be created with this command: 

.. code:: bash

    stemcnv-check setup-files


.. _steps-ext-config:

Setting up the config file
==========================

The command ``stemcnv-check setup-files`` generated a basic config file template containing only the required options. 

The default config file (config.yaml) created by this commands defines all settings required to start the analysis and 
inherits other settings from the inbuilt default. By using the command line option ``--config-details`` an extended 
config template containing more options can be generated. Possible parameters are ``medium``, ``advanced`` or ``complete``. 
**Altering these advanced parameters is not recommend for standard use as they can change analysis outcome!**
| The complete config is described in the :ref:`technical detail section <tech-config>`
  
**Edit the config file** so that all entries marked as ``“#REQUIRED”`` are filled in.
  
The config file (default: config.yaml) defines all settings for the analysis and inherits from the inbuilt default, as
well as system-wide array definitions if those exist. While most of the settings can be left on default, the input files
need to be defined. The file paths for these files need to be entered in the config under the 'array_definition' section.

**Array definition**
  
In this section you also need to give your array a name (that needs to match the 'Array_Name' column in the sample table) and define a genome version (hg19 or hg38). 

.. caution::

    Please note that the Illumina bpm and csv manifest files are also specific to a certain genome version, usually files for hg19 end in 'A1' and those for hg38 end in 'A2' (the egt cluster file is not specific and can be used for any genome version).

Other array specific files mentioned in the config can be auto-generated (see next step below).
While most of the settings can be left on default, the input files need to be defined. Among those are also the files for the definition of the array platform, which are the primary
required settings apart from raw data locations, that can not have defined defaults in the config file created by the
setup-files command.

- **'ExampleArray'** should to be renamed to the actual array name
- **genome_version** options are: hg38/GRCh38 or hg19/GRCh37

**Define files specific to the used array platform and genome build:**

- **egt_cluster_file**: the illumina cluster file (.egt) for the array platform, available from Illumina or the provider running the array

- **bpm_manifest_file**: the beadpool manifest file (.bpm) for the array platform, available from Illumina or the provider running the array

- **csv_manifest_file** (optional): the manifest file in csv format, available from Illumina or the provider running the array

- **raw_data_folder**: input folder, path to the input directory under which the raw data (.idat) can be found. 
  This folder should contain subfolders that match the Chip_Name column in the sample table (containing the array chip IDs).
  **idat files should be grouped in a subfolder per array-chip (sentrix_name).**

- **data_path**: the output of StemCNV-check will be written to this path
- **log_path**:  output folder, stemcnv-check will write log filesthe log files of StemCNV-check to this path

.. code:: yaml

  array_definition:
    GSAMD-24v3-0:
      genome_version: 'hg19'
      bpm_manifest_file: '../cluster-manifest-data/GSAMD-24v3/gh19/GSAMD-24v3-0-EA_20034606_A1.bpm'              
      egt_cluster_file: '../cluster-manifest-data/GSAMD-24v3/gh19/GSAMD_24v3-0_A1-LAB-2235HiQ-Samples.egt'    
      csv_manifest_file: '../cluster-manifest-data/GSAMD-24v3/gh19/GSAMD-24v3-0-EA_20034606_A1.csv'
      penncnv_pfb_file: '__cache-default__'
      penncnv_GCmodel_file: '__cache-default__'
      array_density_file: '__cache-default__'
      array_gaps_file: '__cache-default__'

  raw_data_folder: '../RAW_DATA'
  data_path: data
  log_path: logs

  reports:
    StemCNV-check-report:
      file_type: 'html'


.. _steps-ext-sampletable:

Setting up the sample table
===========================

The sample table (default: sample_table.tsv) is a tab-separated file describing all samples to be analyzed.
**Excel or tsv** formats are supported.

The default format of the sample table is tsv. If you prefer to use an xlsx file here you can create an example by using:

``stemcnv-check setup-files --sampletable-format xlsx``

You can also use your own Excel file, if the following criteria are met:

  - The actual sample table is in the first sheet of the file and this sheet *only* contains columns for the sample table (optionally with commented lines starting with a '#')

  - All required columns are present and correctly named (the order of columns is not important)

  - It is possible to deviate from the standard column names, but the expected column names need be contained in the 
    actual column names and there needs to a singular way to extract them (via regex).

      - In this case you need to use the ``--column-remove-regex`` option to tell the pipeline how to modify your column 
        names to derive the expected names. If used without an explicit regex (for expert users) spaces and anything 
        following them will be removed from your column names.

      - A simple example with ``--column-remove-regex`` (default) option would be to use i.e:
        'Sample_ID for pipeline', 'Chip_Name (Sentrix Barcode)', 'Chip_Pos (Sentrix Position)'

Filling in the sample table with your data
------------------------------------------

| **Required Columns**:
| *Sample_ID, Chip_Name, Chip_Pos, Array_Name, Sex, Reference_Sample, Regions_of_Interest, Sample_Group*
|
| Specific explanations for columns:

  - Sample_ID

    The folder and file names for samples are derived from this entry. All entries *must* be unique. 
    To prevent issues with filenames only alphanumeric characters (all letters and number) and the characters ``-`` 
    and ``_`` (dash and underscore) are allowed by default.

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

    This column can define sample specific regions of interest (i.e. gene edited sites), if none exist it can be left 
    empty but still must be included in the sample sheet.  
    The syntax for regions of interest is ``NAME|region``, the ``NAME|`` part is optional and mainly useful for 
    labeling or describing the region.  
    The ``region`` part is mandatory and can be one of the following:

    1) **Position:** "chrN:start-end": ``chrN`` can be i.e. ``chr3`` or just ``3``, start and end are coordinates (which are genome build specific!)
    2) **Genomic band:** i.e. "4q21.3": a cytogenetic band, both full bands (q21) and subbands (q21.3) are allowed 
    3) **Gene symbol:** i.e. "TP53": The gene name (or symbol) needs to exactly match the reference annotation (UCSC gtf). Validity of gene symbols can be used using the `HGNC Multi-symbol checker <https://www.genenames.org/tools/multi-symbol-checker/>`_
    
    Multiple regions for a single sample should all be in one column entry and be separated by a ``;``

  - Sample_Group

    This column can be used for grouping of related samples. By default all samples within the same group will be included in sample comparison based on SNP clustering.

								
.. list-table::  Example Sample table
   :widths: 15 15 10 10 10 10 10 10
   :header-rows: 1
								
   * - Sample_ID 
     - Chip_Name
     - Chip_Pos
     - Array_Name
     - Sex
     - Reference_Sample
     - Regions_of_Interest
     - Sample_Group
   * - HG001
     - 207521920117
     - R09C02
     - ExampleArray
     - female
     - donor_fibroblasts_HG001
     -
     - Group1
   * - donor_fibroblasts_HG001
     - 207521920117
     - R05C02
     - ExampleArray
     - male
     -
     -
     - Group1
   * - hESC_1
     - 207521920117
     - R07C02
     - ExampleArray
     - female				
     -
     - 4q21.3
     - 
   * - HG005
     - 207521920117
     - R01C02
     - ExampleArray
     - male
     -
     -
     - HG006
   * - HG006
     - 207521920117
     - R03C02
     - ExampleArray
     - male
     -
     -
     - 
   * - HG007
     - 207521920117
     - R11C02
     - ExampleArray
     - female
     -
     -
     - 


.. _steps-ext-staticdata:

Static files generation
=======================

This step takes place after the  sample data for that array is available, sample table and the config file have been set up.

**Array & genome-build specific static files** are automatic generated.

.. code:: bash

   stemcnv-check make-staticdata


.. note::

    This step will also include **download of fasta and gtf** file for the reference genome build.**
    Array specific files and an updated array_definition block for the config will be written into the 
    cache directory (default: `~/.cache/stemcnv-check`).


StemCNV-check generally requires two types of static data files: those that are specific to the genome version (incl. 
the genome reference sequence) and those that are specific to the array platform. All of these files can be downloaded 
or generated by StemCNV-check using the ``stemcnv-check make-staticdata`` command, however array specific files can only 
be created if raw data for at least one sample is available. Usually genome version specific files are only downloaded 
once and saved in a central cache location, so they should already be available after running the example data.  
The files specific to an array platform are also saved to this central cache, so that they can be shared between different 
projects. Additionally, an updated array definition block for the config is written to the cache, so that the array 
definition is also saved. However, array definitions from a project specific config file will still take precedence over 
the central definitions.

To create the array specific files, follow these steps: 

- make sure that the sample table and config file, with all required entries, are correctly set up
- Run the ``stemcnv-check make-staticdata`` 

  - This command will download missing genome specific files from the internet

  - Then it will generate the array specific files, which also requires processing the raw data from at least one sample.

  .. tip:: 
    If you already have a genome reference fasta on your system you can also use that, 
    instead of downloading a second one. To do so you need to provide the path to the fasta file for the corresponding 
    genome version in the 'global_settings' block of the config file. This section will only be included in the config 
    if you use at least the ``--config-details medium`` flag for the setup-files command.

This command will also print out the paths to the generated array specific files. You can either copy these paths your 
project specific config file to use a complete array definition, or you can simply remove the array definition block 
and rely on the automatically saved central definitions.

