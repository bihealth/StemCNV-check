.. _steps3-basic-usage:

Basic Usage
^^^^^^^^^^^

This section explains the basics of how to use StemCNV-check once it is installed *without* going into much detail.  
If you want detailed instructions for running an :ref:`example dataset <steps4-example-data>` or setting up a 
:ref:`new project from scratch <steps-ext-project-setup>`, please take a look our respective step-by-step instructions.

**Setting up a StemCNV-check analysis project requires:**

- array raw data files of samples to be included in analysis
- manifest and cluster files of the specific array used: *egt_cluster_file, bpm_manifest_file, csv_manifest_file (optional)*
- config file
- sample table

It is generally recommended to keep each analysis project and it's associated files in a separate folder.

Before running StemCNV-check, you will need to create and fill out a **sample-table** and **config file** in the 
project folder and then create static data for the specific SNP array used.  

Create a project folder using the following command

.. code-block:: bash

    mkdir {project_name}


Setting up the sample-table and config files
============================================

On the comand line change into the project folder:

.. code-block:: bash

    cd {project_name}

Empty example files for the *sample table* and *config file* can then be created with this command:

.. code-block:: bash

    stemcnv-check setup-files

If you prefer to use an xlsx (Excel) file for the sample table it can be created by using:

.. code-block:: bash

    stemcnv-check setup-files --sampletable-format xlsx


.. tip::

    StemCNV-check will by default always look for a "sample_table.tsv" (or .xlsx) and a "config.yaml" file.
    However, files with different names can also be used (i.e. to compare the results of different settings), 
    but then need to be passed via the ``--sample-table`` or ``--config`` command line options, 
    which work for all StemCNV-check commands.


Sample Table
------------

The **sample table** describes all samples you want to analyse in a given project. 
It can be expanded later in a project to add more samples, without affecting the primary analysis of older samples 
(though some comparisons may be updated). The following columns of the sample table are required for running StemCNV-check:

*Sample_ID, Chip_Name, Chip_Pos, Array_Name, Sex, Reference_Sample, Regions_of_Interest, Sample_Group*

The first 5 columns (Sample_ID - Sex) are required to be filled out for each sample. 
The columns *Reference_Sample* (used to define the sample to compare with i.e. originating fibroblast or parental cell line) 
and *Sample_Group* (used for identity comparison) are optional and have be left blank, if no fitting values can be entered.

The sample table file created by the ``setup-files`` command contains comments (lines starting with a hash ``#`` symbol, 
which are ignored by StemCNV-check), explaining the columns in more details. These can be removed from the sample table.

For more in-depth instructions and explanations about the sample table see :ref:`Sample Table setup details <steps-ext-sampletable>`.

Config File
-----------
The **config file** contains all settings of StemCNV-check. By default, the config file created by the ``setup-files`` 
command only has the minimum number of entries that are required for for running an analysis project. These entries that 
need to be filled in are marked with a ``#REQUIRED”`` comment. All other (optional) settings are instead taken from inbuilt defaults. 

These include specifically the file paths to array manifest files (describing the array probes) and the input and output 
file paths the pipeline should use:

- **'ExampleArray'** should to be renamed to the actual array name
- **genome_version:** can be set to hg38/GRCh38 or hg19/GRCh37. 

Please note that the Illumina bpm and csv manifest files are also specific to a certain genome version, usually files for hg19 end in ‘A1’ and those for hg38 end in ‘A2’ (the egt cluster file is not specific and can be used for any genome version)

- **egt_cluster_file**: the illumina cluster file (.egt) for the array platform, available from Illumina or the provider running the array
- **bpm_manifest_file**: the beadpool manifest file (.bpm) for the array platform, available from Illumina or the provider running the array
- **csv_manifest_file** (optional): the manifest file in csv format, available from Illumina or the provider running the array
- **raw_data_folder**: input folder, path to the input directory under which the raw data (.idat) can be found. Ths folder should contain subfolders that match the Chip_Name column in the sample table (containing the array chip IDs). **idat files should be grouped in a subfolder per array-chip (sentrix_name).**

An example configuration that assumes the sub folders "cluster-manifest-data" (containg the cluster and manifest files) and "raw_data" (containing the array raw data) is given below:

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
    
    raw_data_folder: '../raw_data' 
    data_path: data
    log_path: logs
    
    reports:
      StemCNV-check-report:
        file_type: 'html'

For more in-depth instructions and explanations about the config file see :ref:`config file setup details <steps-ext-config>` .

Generating array static data
----------------------------

StemCNV-check requires some array specific additional files that are separate from the array manifests, but are also 
*static*, i.e. they only need be created once. Some of these files require information that is only accesible after 
pre-processing at least one sample, so you need a filled out config file and sample table first.
StemCNV-check has an inbuilt workflow to create these files that also saves these files independently from the 
project, so they can be re-used later (see :ref:`file caching <tech-cache>`). This requires that the same ``ArrayName`` 
is used in the sample table (and config) file across different projects.
In addition, the same workflow will also download other information required, like the genome reference files.

The workflow to create all static files and prepare StemCNV-check can be started with this command:

.. code-block:: bash

    stemcnv-check make-staticdata


.. tip::

    If you also run other bioinformatics analysis, you may already have genome ``fasta`` and ``gtf`` files on your system.
    In this case, you can configure StemCNV-check to use those files instead of downloading new ones. This needs to be set
    in the ``global_settings`` part of the config, which is included from ``--config-details medium`` and above.

Starting the StemCNV-check analysis
-----------------------------------

After config file and sample-table have been set up and the static data for an array has been created, the StemCNV-check 
workflow can be started with this command:

.. code-block:: bash

    stemcnv-check run


.. tip::

    StemCNV-check is built on snakemake and can also utilise all of snakemake's advanced features. 
    You can forward command line options to snakemake by separating them with a ``--``. This way you can for example 
    make use of snakemake executors that can interface with HPC scheduling systems: ``stemcnv-check run -- --executor slurm``

| After the analysis finised successfully reports can be found in the foder defind in the config file ``data_path`` 
| eg.: ``./data/{sample name}/{sample name}.StemCNV-check-report.html``

