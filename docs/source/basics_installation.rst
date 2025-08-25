.. _basics-installation:

Installation
^^^^^^^^^^^^

Ths section contains step-by-step instructions for setting up StemCNV-check. All commands should be run in a linux 
terminal, including WSL (Windows Subsystem for Linux) on Windows.

.. tip::

    If you are not familiar with using scientific software on Linux and/or are a Windows user, it is strongly 
    recommended to read through the :ref:`Basics for non-linux Users <basics-nonlinux>` section first.  
    This will also explain how to setup WSL (Windows Subsystem for Linux) on Windows.

Installation of Conda
======================

.. note:: 

    Conda is a software that facilitates the distribution and installation of primarily scientific software with the 
    ability to control which specific versions are used. StemCNV-check utilises this for almost all steps of the 
    workflow and as such depends on a working conda setup. In principle any conda setup can be used, but for anyone 
    not familiar we recommend the following.

Install the `miniforge conda <https://github.com/conda-forge/miniforge>`_.  
In short, use these commands: 

.. code:: bash

    wget "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
    bash Miniforge3-$(uname)-$(uname -m).sh

| During the installation follow instructions and suggestions that are displayed in the terminal. When the istaller askes if it should update the shell profile to the allow automatic initialisation of conda, answer/type: `yes` instead of the default no. In all other cases the default answers should let you continue.
| After successful installation of conda WSL has to be restarted (close and reopen the window) or reload the updated 
.bashrc of WSL using the command: ``source ~/.bashrc``


Installation of StemCNV-check
=============================

It is recommended to install StemCNV-check with conda through the bioconda channel. 

Install the latest released version through conda (**recommended**)
------------------------------------------

To install `StemCNV-check from bioconda <https://anaconda.org/bioconda/stemcnv-check>`_, run this command:

.. code:: bash

    conda install bioconda::stemcnv-check
   


Install the development version from github (**only for developers**)
-------------------------------------------
.. note:: 

    If you are an experienced conda user, you can of course create a specific environment for StemCNV-check.

Alternatively to bioconda, one can perform the following steps to install and setup up the development version of StemCNV-check.

- Clone the StemCNV-check git repository

.. code:: bash

   git clone https://github.com/bihealth/StemCNV-check.git

- Set up conda environment for StemCNV-check and install dependencies

  - for WSL (on Windows)

    .. code-block:: bash

      conda install python=3.12 "gcc_linux-64<14" apptainer fuse-overlayfs

  - for Linux

    .. code-block:: bash

      conda install python=3.12

.. tip::
  If you also use conda for other projects, you may prefer to use a specific environment only for StemCNV-check:
  ``conda create -n stemcnv-check python=3.12; conda activate stemcnv-check``


- Change into the StemCNV-check directory:

  .. code-block:: bash

    cd StemCNV-check

- Install StemCNV-check and its dependencies with pip:

  .. code-block:: bash

    pip install -e .

**Updating the developement version**

As long as you are in the StemCNV-check directory you can update the development version of StemCNV-check with this 
command:

.. code-block:: bash

  git pull
  pip install -e .
