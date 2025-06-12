
Setting up the config file
============

The default config file (config.yaml) defines all settings for the analysis and inherits from the inbuilt default.
  
**Edit the config file** so that all entries marked as
``“#REQUIRED”`` are filled in.
  
The config file (default: config.yaml) defines all settings for the analysis and inherits from the inbuilt default, as
well as system-wide array definitions if those exist. While most of the settings can be left on default, the input files
need to be defined. The file paths for these files need to be entered in the config under the 'array_definition' section.

**Array definition**
  
In this section you also need to give your array a name (that needs to match the 'Array_Name' column in the sample table) and define a
genome version (hg19 or hg38). Please note that the Illumina bpm and csv manifest files are also specific to a certain
genome version, usually files for hg19 end in 'A1' and those for hg38 end in 'A2' (the egt cluster file is not specific
and can be used for any genome version).
Other array specific files mentioned in the config can be auto-generated (see next step below).
While most of the settings can be left on default, the input files need to be defined. Among those are also the files for the definition of the array platform, which are the primary
required settings apart from raw data locations, that can not have defined defaults in the config file created by the
setup-files command.

- **'ExampleArray'** should to be renamed to the actual array name

- **genome_version options:** hg38/GRCh38 or hg19/GRCh37
**Define  files specific to the used array platform and genome build:**

- **egt_cluster_file**: the illumina cluster file (.egt) for the array platform, available from Illumina or the provider running the array

- **bpm_manifest_file**: the beadpool manifest file (.bpm) for the array platform, available from Illumina or the provider running the array
- **csv_manifest_file** (optional): the manifest file in csv format, available from Illumina or the provider running the array

- **raw_data_folder**: input folder, path to the input directory under which the raw data (.idat) can be found. Ths folder should contain subfolders that match the Chip_Name column in the sample table (containing the array chip IDs). **idat files should be grouped in a subfolder per array-chip (sentrix_name).**

- **data_path**: the output of StemCNV-check will be written to this path
- **log_path**:  output folder, stemcnv-check will write log filesthe log files of StemCNV-check to this path

.. code:: bash

   array_definition:
        GSAMD-24v3-0:
          genome_version: 'hg19'
          bpm_manifest_file: 'static-data/ExampleArray/GSAMD-24v3-0-EA_20034606_A1.bpm'
          csv_manifest_file: 'static-data/ExampleArray/GSAMD-24v3-0-EA_20034606_A1.csv.gz'
          egt_cluster_file: 'static-data/ExampleArray/GSAMD-24v3-0-EA_20034606_A1.csv.gz'
          #Optional (leave empty if not used)
          penncnv_GCmodel_file: 'static-data/ExampleArray/GSAMD-24v3-0-EA_20034606_A1.egt'
          array_density_file: 'static-data/ExampleArray/density_hg19_ExampleArray.bed'
          array_gaps_file: 'static-data/ExampleArray/gaps_hg19_ExampleArray.bed'
          penncnv_pfb_file: 'static-data/ExampleArray/PennCNV-PFB_hg19_ExampleArray.pfb'

   raw_data_folder: ../RAW_DATA
   data_path: data_scoring
   log_path: logs

   reports:
        StemCNV-check-report:
          file_type: 'html'



**Advanced options**

The config file created by this command will only include the absolute necessary settings to run the workflow. If
you are interested in setting additional parameters or changing the content of the report, you can add this flag
--config-details medium to the command (also available with ‘advanced’ or ‘complete’ instead of ‘medium’)
