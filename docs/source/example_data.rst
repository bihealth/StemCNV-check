
Running the example dataset
=============
  
**Example data** (using data from the Genome in a Bottle samples) can be used to test the setup. 
  
Test data can be downloaded via git LFS and StemCNV-check can be run with the following commands. 
  
Obtain the example data:
-------------------
Note that this will download a fasta and gtf file for the human genome. If you have suitable files available locally, it
is recommended to replace the corresponding paths in the config.yaml to avoid unnecessary and time-consuming
downloads.
  
• Install git-lfs: ``sudo apt-get update; sudo apt-get install git-lfs``
  
• Download the test data:
  
  – go to the StemCNV-check directory, if you are not already there: ``cd StemCNV-check``
  
  – Actually download the files (~300Mb): ``git lfs fetch; git lfs checkout``
  
Run the example data:
-------------------

• First change into the example_data folder: ``cd example_data``
  
• Now prepare additional static-data:
``stemcnv-check make-staticdata``
  
– This will download and generate several files, including the human reference genome, so it requires some
disk space and a bit of time.
  
– Note: experienced users can also provide paths to the necessary fasta and gtf files from their own system.
This requires entries in the global_settings section of the config.yaml file, which is included with medium
detail when creating example files.
  
• **Finally start the StemCNV-check workflow for the example data:** ``stemcnv-check``
  
– **Important for WSL users:** stemcnv-check (or rather the underlying snakemake process) will by default
attempt to use all available memory. On WSL this can easily lead to cases where windows (or
WSL) start stopping indivdual processes due to memory usage. To prevent this, limit the overall
memory usage by instead using the following command to start the workflow. You can adjust the
memory limit to your system, this assumes that at least 6Gb of memory are available for the workflow:
``stemcnv-check run -- --resources mem_mb=6000``

– Note in case of (download/network) issues: please see section 4.4 for potential solutions to common issues.
