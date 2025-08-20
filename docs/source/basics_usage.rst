
Basic Usage
^^^^^^^^^^^

Essentially: 

::include:: ../../README.md

First analysis
============
Before the first analysis sample table and config file need to be set up (see above). Unless otherwise specified, stemcnv-check defaults to look for a "sample_table.tsv" (or .xlsx) and "config.yaml" file.

It is recommended to start by **creating a separate folder** for your project. This folder should include raw data folder, config.yaml and sample table files.

Config file settings
============

Setting analysis parameters or changing them requires editing the text in the generated default config file. Start by opening the config.yaml in text editor. Then type in or change the necessary parameters. 

The default config file (config.yaml) defines all settings for the analysis and inherits from the inbuilt default.

Adjust the config file so that all entries marked as ``“#REQUIRED”`` are filled in.


Continuing from previous steps and running multiple projects
------------------------------------------------

If you have already completed some or all steps of the initial installation setup, you do not need to repeat them. 
However, if you restart the terminal (WSL) window, this will also always return you to the home directory. 
To make sure that you are in the correct directory for the following steps, each section assumes you start (again) in
the home directory. You can restart your WSL window or use cd ~ to get back to the home directory. 
If you follow the recommended instructions the stemcnv-check command will always be available in your terminal
(through the base environment of conda), but if you use the version for experienced users you will need to activate
your stemcnv-check conda every time you start a new terminal window: conda activate stemcnv-check.

If you want to run your own data analysis project(s) with StemCNV-check, it is recommended to create a new
directory for each project. These directories can be located anywhere on your system, you can create them in the
home directory (of WSL for windows) with the following command: mkdir {project_name} . You can also use
any directory/folder in your normal Windows/MacOS system, however in that case you will need to either open
the (WSL) terminal in that directory or change to if normally starting WSL (using cd ). On windows open WSL
in a specific folder can be done by right-clicking in the folder and selecting ‘Open in Terminal’ or ‘Open in WSL’
(potentially listed under ‘Show more Options’).