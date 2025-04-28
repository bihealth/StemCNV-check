Running StemCNV-check and project setup
============


**Setting up your own data for analysis with StemCNV-check requires:**

- config file
- manifest and static files: egt_cluster_file, bpm_manifest_file, csv_manifest_file (optional)
- sample table



Config file
-----------
The config file (default: config.yaml) defines all settings for the analysis and inherits from the inbuilt default, as 
well as system-wide array definitions if those exist. While most of the settings can be left on default, the input files 
need to be defined. Among those are also the files for the definition of the array platform, which are the primary 
required settings apart from raw data locations, that can not have defined defaults in the config file created by the 
setup-files command. The array definition needed for StemCNV-check includes the following files:

The default config file (config.yaml) defines all settings for the analysis and inherits from the inbuilt default.

*Adjust the config file* so that all entries marked as
``“#REQUIRED”`` are filled in.

**Define  files specific to the used array platform and genome build:**

- **egt_cluster_file**: the illumina cluster file (.egt) for the array platform, available from Illumina or the provider running the array

- **bpm_manifest_file**: the beadpool manifest file (.bpm) for the array platform, available from Illumina or the provider running the array
- **csv_manifest_file** (optional): the manifest file in csv format, available from Illumina or the provider running the array


The file paths for these files need to be entered in the config under the 'array_definition' section. In this section 
you also need to give your array a name (that needs to match the 'Array_Name' column in the sample table) and define a 
genome version (hg19 or hg38). Please note that the Illumina bpm and csv manifest files are also specific to a certain 
genome version, usually files for hg19 end in 'A1' and those for hg38 end in 'A2' (the egt cluster file is not specific 
and can be used for any genome version).  
Other array specific files mentioned in the config can be auto-generated (see next step below).

**The config file needs to define the following paths:**

- **raw_data_folder**: path to the input directory under which the raw data (.idat) can be found. Ths folder should contain subfolders that match the Chip_Name column in the sample table (containing the array chip IDs)

- **data_path**: the output of StemCNV-check will be written to this path
- **log_path**: the log files of StemCNV-check will be written to this path

**Advanced options**

The config file created by this command will only include the absolute necessary settings to run the workflow. If
you are interested in setting additional parameters or changing the content of the report, you can add this flag
--config-details medium to the command (also available with ‘advanced’ or ‘complete’ instead of ‘medium’)


Sample table
-----------

The sample table (default: sample_table.tsv) is a tab-separated file describing all samples to be analyzed. Excel or tsv formats are supported.
Empty example files for the sample table and config can be created with this command:

``stemcnv-check setup-files``

If you prefer to use an xlsx file here you can create an example by using:  

``stemcnv-check setup-files --sampletable-format xlsx``

**Fill in the sample table with your data**

- Required columns: Sample_ID, Chip_Name, Chip_Pos, Array_Name, Sex, Reference_Sample
- Optional (reserved) columns: Regions_of_Interest

You can also use your own Excel file, if the following criteria are met:

- The actual sampletable is in the first sheet of the file and this sheet *only* contains columns for the sample table 
  (optionally with commented lines starting with a '#')

- All required columns are present and correctly named (the order of columns is not important)
  - It is possible to deviate from the standard column names, but the expected column names need be contained in the acutal 
    column names and there needs to a singular way to extract them (via regex). 
  - In this case you need to use the ``--column-remove-regex`` option to tell the pipeline how to modify your column names 
    to derive the expected names. If used without an explicit regex (for expert users) spaces and anything following 
    them will be removed from your column names.

  - A simple example with ``--column-remove-regex`` (default) option would be to use i.e:  
    'Sample_ID for pipeline', 'Chip_Name (Sentrix Barcode)', 'Chip_Pos (Sentrix Position)'

								
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

**Extended sample table. Description of the data types contained in the columns.**

- Sample_ID
Include bank ID when possible, only: - or _, do not use special characters: (), {}, /, \, ~,*, & Name has to be UNIQUE.
This column has auto-formatting enabled, so that the IDs will work with the CNV-pipeline:
	- red entries are either duplicate or contain not-allowed characters (/ and .\)
	- orange entries contain characters that the pipeline will remove (since they can cause issues if used in file names):  :,;()[]{}!?* and <space>

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

Generating (array specific) static files
-----------------------------------

StemCNV-check generally requires two types of static data files: those that are specific to the genome version (incl. 
the genome reference sequence) and those that are specific to the array platform. All of these files can be downloaded 
or generated by StemCNV-check using the ``stemcnv-check make-staticdata`` command, however array specific files can only 
be created if raw data for at least one sample is available. Usually genome version specific files are only downloaded 
once and saved in a central cache location, so they should already be available after running the example data.  
The files specific to an array platform are also saved to this central cache, so that they can be shared between different 
projects. Additionally, an updated array definition block for the config is written to the cache, so that the array 
definition is also saved. However, array definitions from a project specific config file will still take precedence over 
the central definitions, therefore the project specific config file need to be adapted once more after generating the
array specific static files.

To create the array specific files, follow these steps: 

- make sure that the sample table and config file, with all required entries, are correctly set up
- Run the ``stemcnv-check make-staticdata`` 
  - This command will download missing genome specific files from the internet
    - if you already have a genome reference fasta on your system you can also use that, 
      instead of downloading a second one. To do so you need to provide the path to the fasta file for the corresponding 
      genome version in the 'global_settings' block of the config file. This section will only be included in the config 
      if you use at least the ``--config-details medium`` flag for the setup-files command. Other files like gtf can also
  - Then it will generate the array specific files, which also requires processing the raw data from at least one sample.

This command will also print out the paths to the generated array specific files. You can either copy these paths your 
project specific config file to use a complete array definition, or you can simply remove the array definition block 
and rely on the automatically saved central definitions.

Running the analysis
-----------

Once the config file is properly set up and the necessary static files are generated, you can run the StemCNV-check
workflow with simple command:   
``stemcnv-check``

*Reminder for WSL users:* as before you may need to limit the memory usage of the workflow 
and use this command or a variantion instead: ``stemcnv-check run -- --resources mem_mb=6000``
