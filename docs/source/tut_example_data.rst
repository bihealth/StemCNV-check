
Running the example dataset
^^^^^^^^^^^^^^^^^^^^^^^^^^^
  
We provide an example data set, based on samples form Genome in a Bottle samples, that can be used to test the setup.
  
Obtain the example data
=======================

.. important::
    The test data is currently only available via github. Therefore this instruction set assumes you have 
    cloned the StemCNV-check repository from github (to do so, follow the installation instructions for the 
    github/developement verion until the ``git clone ...`` command).  
    We will make the test data available as a data resource in the future.

Test data can be downloaded via git LFS and StemCNV-check can be run with the following commands. 


• Install git-lfs: ``sudo apt-get update; sudo apt-get install git-lfs``
  
• Download the test data:
  
  – go to the StemCNV-check directory, if you are not already there: ``cd StemCNV-check``
  
  – Actually download the files (~300Mb): ``git lfs fetch; git lfs checkout``
  
Run the example data
====================

• First change into the example_data folder: ``cd example_data``
  
• Now prepare additional static-data: ``stemcnv-check make-staticdata``
  
  .. note:: 
    This will download a fasta and gtf file for the human genome. This download will require some disk space and 
    a bit of time.  
    If you are an experienced user, you can also provide paths to the necessary fasta and gtf files from your own system.
    This requires entries in the global_settings section of the config.yaml file, which is included with medium
    detail when creating example files.

  .. tip::
    If StemCNV-check encouters issues automatically download files, please see the Troubleshooting section 
    for potential solutions to common issues.

• Finally start the StemCNV-check workflow for the example data: ``stemcnv-check run``

