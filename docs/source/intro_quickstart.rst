.. _intro-quickstart:

Quickstart for experienced users
================================

Installation
------------

StemCNV-check can be installed via `bioconda <https://bioconda.github.io/recipes/stemcnv-check/README.html>`_ or 
from the development version on `github <https://github.com/bihealth/StemCNV-check>`_:

.. code:: bash

    conda install bioconda::stemcnv-check

StemCNV-check requires at least 8GB of RAM and 4 CPU cores with default settings, however more is recommended and 
parallelization through snakemake can efficiently utilize more resources form more powerful machines.

Setup
-----

Before running StemCNV-check, you will need to create and fill out a **sample-table** and **config file**.
Example files for the *sample table* and *config file* with comment annotations can be created with this command:

.. code-block:: bash

    stemcnv-check setup-files [--config-details {level}]

StemCNV-check needs to create some static data for each SNP array. Apart from the manufacturer provided array manifest files,
these can be automatically generated. This requires a filled out config and at least one sample:

.. code-block:: bash

    stemcnv-check make-staticdata

.. tip::

    This command will also download fasta and gtf files for the select genome version. If you already have some on your system, 
    you can also configure StemCNV-check to use those files instead of downloading new ones. This needs to be set
    in the ``global_settings`` part of the config, which is included from ``--config-details medium`` and above.

    Alternatively you can create links pointing to matching files under these paths (see also :ref:`file caching <tech-cache>`):

    - Tabix indexed fasta: ``~/.cache/stemcnv-check/fasta/homo_sapiens/113_GRCh{37,38}/Homo_sapiens.GRCh{37,38}.dna.primary_assembly.fa.gz``
    - gtf: ``~/.cache/stemcnv-check/static-data/gencode.hg{19,38}.v45.gtf.gz``
      

More in-depth instructions for these setup steps can be found in the :ref:`extended project setup <steps-ext-project-setup>` section.

Analysis
--------

Once the setup for a specific array is done, the analysis workflow can be started with a single command:

.. code-block:: bash

    stemcnv-check run


.. tip::

    StemCNV-check is built on snakemake and can also utilise all of snakemake's advanced features. 
    You can forward command line options to snakemake by separating them with a ``--``. This way you can for example 
    make use of snakemake executors that can interface with HPC scheduling systems: ``stemcnv-check run -- --executor slurm``


New samples can simply be added to the sampletable and will be analysed with the next run command. 
Similarly, after changing options in the config files the workflow will update and rerun all affected steps automatically.

A full list of all config options can be found in the :ref:`config file reference <tech-config>`.