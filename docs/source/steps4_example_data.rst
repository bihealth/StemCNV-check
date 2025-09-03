.. _steps4-example-data:

Running the example dataset
^^^^^^^^^^^^^^^^^^^^^^^^^^^
  
We provide an `example data set <https://zenodo.org/records/16962381>`_, based on samples form Genome in a Bottle samples, 
that can be used to test the setup.

Please note that the Genome in a Bottle samples are *not* pluripotent stem cells, so interpretation of the data is not really possible.
  
Obtain the example data
=======================

The example input data can be downloaded from `zenodo <https://zenodo.org/records/16962381/files/example_input_data.zip?download=1>`_,
or directly on the command line:
  
  .. code-block:: bash
  
    wget https://zenodo.org/records/16962381/files/example_input_data.zip
    unzip example_input_data.zip

If you download the data through the web-browser make sure to unzip it also. 
Even when using WSl, it is recommended to leave all StemCNV-check files on the linux file system to improve performance.

  
Run the example data
====================

.. note::

    Working with the example data requires, that you have :ref:`installed StemCNV-check <steps2-installation>`.

The example data fileset contains all necessary input files as well as a fully filled in sample table and config file, 
so you can proceed to generation of static data and data analysis right away:


- Make sure you have an open terminal, located in unpacked the example data folder containing the config.yaml and sample table files.
  If you used the above commands to download the files from the command line this should already be the case.

- Start the preparation of the static data for the example data

  .. code-block:: bash

    stemcnv-check make-staticdata
  
  .. note:: 
    If you have not used StemCNV-check before, this step will include the download of reference fasta and gtf files for 
    the human genome. This download will require some disk space and a bit of time.  

    If you are an experienced user, you can also provide paths to the necessary fasta and gtf files from your own system.
    This requires entries in the global_settings section of the config.yaml file, which is included with medium
    detail when creating example files.

  .. tip::
    If StemCNV-check encounters issues automatically download files, please see the Troubleshooting section 
    for potential solutions to common issues.

- Finally start the StemCNV-check workflow for the example data: 

  .. code-block:: bash

    stemcnv-check run

After performing all these steps you should have complete output for all 6 samples in the ``data`` folder.

You can also download the `example output data <https://zenodo.org/records/16962381/files/example_output_data_v1.0.zip?download=1>`_, 
created with StemCNV-check version 1.0, to compare them with your newly generated files.
