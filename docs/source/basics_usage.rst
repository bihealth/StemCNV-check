.. _basics-usage:

Basic Usage
^^^^^^^^^^^

.. admonition:: Under construction

    This page is still under construction and has not been finalised yet

This section explains the basics of how to use StemCNV-check once it is installed *without* going into much detail.  
If you want detailed instructions for running an :ref:`example dataset <tut-example-data>` or setting up a 
:ref:`new project from scratch <tut-project-setup>`, please take a look our tutorials.


Before running StemCNV-check in any project, you will need to setup and fill out the sample table and config file 
specific for that project and then create rge array specific static data.  
It is generally recommended to keep each project and it's associated files in a separate folder. 
Additionally, StemCNV-check will by default always look for a "sample_table.tsv" (or .xlsx) and a "config.yaml" file.
However, different files can also be used (i.e. to compare the results of different settings), but then need to be 
passed via the ``--sample-table`` or ``--config`` command line options, which work for all StemCNV-check commands.

Setting up the sample-sheet and config files
--------------------------------------------

You can create (empty and minimalistic) examples for these files with the following command:

.. code-block:: bash

    stemcnv-check setup-files


The **sample table** describes all samples you want to analyse in a given project. It can be expanded later in a project 
to add more samples, without affecting the primary analysis of older samples (though some comparisons may be updated).  
The following columns of the sample table are required for running StemCNV-check:
- Sample_ID, Chip_Name, Chip_Pos, Array_Name, Sex, Reference_Sample, Regions_of_Interest, Sample_Group
The first 5 of these (Sample_ID - Sex) are required for all samples, Reference_Sample is used to track the origin of a 
sample (i.e. originating fibroblast or master bank) and should be used where possible, the last two columns can be filled optionally.
The sample table file created by the ``setup-files`` command contains comments (lines starting with a hash ``#`` symbol,
which are ignored by StemCNV-check), explaining the columns in more details. Our :ref:`project setup tutorial <tut-project-setup>` 
also contains more in-depth instructions and explanations about the sample table.

The **config file** contain all settings for StemCNV-check. By default, the config file created by the ``setup-files`` 
command only has the minimum number of entries that are required or recommended for analysis. All other settings are 
instead taken from inbuilt defaults. However, you can use the ``--config-details`` option to include more settings.  
All entries in the config file, that need to be filled before StemCNV-check can be run are marked with a ``#REQUIRED”`` comment. 
These include specifically the file paths to array manifest files (describing the array probes) and the input and output 
file paths the pipeline should use:

- array_definition / ``ArrayName`` / genome_version

  StemCNV-check can operate on both hg38 and hg19. The provided array manifest files (specifically ``.bpm`` and ``.csv``) need to match this. 

- array_definition / ``ArrayName`` / bpm_manifest_file, csv_manifest_file, egt_cluster_file

  The illumina cluster file (.egt) for the array platform, available from Illumina or the provider running the array. 
  The array manifest files, usually available from Illumina or the laboratory performing the array analysis.

- raw_data_folder

  The path to the input directory under which the raw data file (.idat) can be found. 
  This folder should contain subfolders that match the Chip_Name column in the sample table (containing the array chip IDs)

- data_path, log_path
  
  The output folders where StemCNV-check will generate its output. By default these will be two folders (``data`` and 
  ``logs``) in the same directory where StemCNV-check is executed.


Generating array static data
----------------------------

StemCNV-check requires some array specific additional files that are separate from the array manifests, but are also 
*static*, i.e. they only need be created once. Some of these files require information that is only accesible after 
pre-processing at least one sample, so you need a filled out config and sampletable first.
StemCNV-check has an inbuilt workflow to create these files that also saves these files independently from the 
project, so they can be re-used later (see :ref:`file caching <tech-cache>`). This requires that the same ``ArrayName`` 
is used in the sample table (and config) file across different projects.
In addition, the same workflow will also download other required, like the genome reference files.

The workflow to create all static files and prepare StemCNV-check can be started with this command:

.. code-block:: bash

    stemcnv-check make-staticdata


.. tip::

    If you also run other bioinformatics analysis, you may already have genome ``fasta`` and ``gtf`` files on your system.
    In this case, you can configure StemCNV-check to use those file instead of downloading new ones. This needs to be set
    in the ``global_settings`` part of the config, which is included from ``--config-details medium`` and above.


Starting StemCNV-check analysis
-------------------------------

After config and sample-table have been set up and the static data for an array has been created, the StemCNV-check 
workflow can be started with this commandÖ

.. code-block:: bash

    stemcnv-check run

.. tip::

    StemCNV-check is built on snakemake and can also utilise all of snakemake's advanced features. 
    You can forward command like options to snakemake by separating them with a ``--``. This way you can for example 
    make use of snakemake executors that can interface with HPC scheduling systems: ``stemcnv-check run -- --executor slurm``
