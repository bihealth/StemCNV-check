Introduction
=============

StemCNV-check is a tool written to simplify copy number variation (CNV) analysis of SNP array data, specifically for quality control of (pluripotent) stem cell lines. StemCNV-check uses snakemake to run the complete analysis from raw data (.idat) up report generation for all defined samples with a single command. Samples need to be defined in a (tabular) sample table and the workflow settings are defined through a yaml file.

StemCNV-check is developed and primarily tested on Linux. Nonetheless, we aim to make it available for users on other operating systems as well. These instructions are primarily intended for users inexperienced with the command line usage and include step-by-step instructions for setting up and running
StemCNV-check. Specifically, it contains instructions for using StemCN-check on Windows, where first an installation of WSL (Windows Subsystem for Linux) needs to be done. StemCNV-check has been tested to run on WSL using these instructions and a standard PC or laptop.

StemCNV-check is not officially supported on MacOS. However, it may be possible to run StemCNV-check on
MacOS, if the given macOS system supports the bioconda channel (which should be the case for x86_64 and ARM64
processors, but NOT for newer hardware with the M1 processors). If you do try this, most installation commands
(except conda) will have to be done with (home)brew instead of the WSL/linux specific apt-get commands.

**How to read and use these instructions:**

• Code blocks represent commands that should be executed in the terminal as is

  – To paste into the Terminal use ‘Ctrl+Shift+V’ - or in WSL a right click - since ‘Ctrl+V’ does not work
  in a linux/WSL terminal.
  
  – In some cases you will need to modify commands or values. This is indicated by red color highlighting
  and curly brackets. You need to modify the text and remove the brackets before finalising the command.
  In cases only certain options are allowed, these are then shown in the format {option1|option2}.
  Example code with {specific option to edit}, sometimes you can choose {option1|option2}

• Notes: indicate important (caveat) information that you should read before executing any associated commands.
| In many cases these will be specific to certain circumstances like your operating system.

• When running your own data you will also need to edit the sample table and config file.

  – In principle, you can use any program to edit these files. Inside the terminal you can use, i.e. the ‘nano’
  program for small edits: nano path/to/filename.txt
  
  – Files in WSL can be opened and edited from the Windows file explorer:
  Searching under ‘Linux’ / ‘Ubuntu’ / ‘home’ / {yourUsername}
  
  – The sample_table(.tsv) can be edited with a text editor or with programs like Excel. It is also possible
  to use an Excel file as sample table, if the first sheet contains the sample table and only the sample table.
  When using Excel be careful of auto-formatting changing any values.

Should you encounter any issues or have questions, please consult the section on troubleshooting and reporting
issues. Especially in cases of issues encountered during the initial workflow installation and setup also consult section regarding common issues and potential solutions.
