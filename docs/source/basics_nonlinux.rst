.. _basics-nonlinux:

Basics for non-linux Users
^^^^^^^^^^^^^^^^^^^^^^^^^^

StemCNV-check is developed and primarily tested on Linux. Nonetheless, we aim to make it available for users on other 
operating systems as well. These instructions are primarily intended for users inexperienced with the command line 
usage and include step-by-step instructions for setting up and running StemCNV-check. 
Specifically, it contains instructions for using StemCN-check on Windows, where first an installation of WSL 
(Windows Subsystem for Linux) needs to be done. StemCNV-check has been tested to run on WSL using these instructions 
and a standard PC or laptop.

StemCNV-check is not officially supported on MacOS. However, it may be possible to run StemCNV-check on
MacOS, if the given macOS system supports the bioconda channel (which should be the case for x86_64 and ARM64
processors, but NOT for newer hardware with the M1 processors). If you do try this, most installation commands
(except conda) will have to be done with (home)brew instead of the WSL/linux specific apt-get commands.

Code examples and text highlights 
=================================

This manual contains code snippets for copying and executing directly. They look like this:

.. code-block:: bash

    echo "Hello World"

.. tip:: To paste commands into a linux terminal use ‘Ctrl+Shift+V’ - or in WSL a right click - since ‘Ctrl+V’ does not 
  work in a linux/WSL terminal.

In some cases there are also ``code examples`` in the text. The same formatting may be used examples of values to enter into files.

.. tip:: Some examples code block may also contain placeholders or optional parameters, this is inidcated by square and 
  curly brackets respectively:

  - ``[optional]`` indicates an optional parameter that can be omitted.
  - ``{placeholder}`` indicates a placeholder that should be replaced with your own value.


Accessing and editing files
===========================

When running your own data you will also need to edit the sample table and config file.
The sample table can either be a tabular text file (.tsv) or an Excel file (.xlsx), the config file is a text based yaml file (.yaml).
In principle, you can use any program to edit these files. 

- For quick edits inside the terminal you can use i.e. the ‘nano’ program.

  .. code-block:: bash
  
    nano config.yaml

- Files in WSL can be opened and edited from the Windows file explorer:

  Searching under ‘Linux’ / ‘Ubuntu’ / ‘home’ / {yourUsername}

- To edit the text based files you can use the windows notepad or any other text editor. 

  More advanced editors like Notepad++ will also highlight the syntax of the config file, which makes it easier to read and edit.

- When using an excel file as sample table be careful of auto-formatting changing any values. 

  Also, the first sheet needs to contain the sample table and nothing else. All filled will be read by StemCNV-check.


Should you encounter any issues or have questions, please first consult the section on :ref:`troubleshooting <issues-common>` 
and :ref:`reporting issues <issues-support>`.


Installation of WSL (Windows Subsystem for Linux)
=================================================

Please consult the `official instructions <https://learn.microsoft.com/en-us/windows/wsl/install>`_ for installing WSL, 
especially if you encounter any problems.

In short:

- Open PowerShell or Windows Command Prompt in administrator mode by right-clicking and selecting "Run as administrator" 
- Enter

  .. code-block:: bash

    wsl --install


- Follow the installation instructions
- You will (likely) be asked to set a username and password for the linux environment. Do remember those.
 
You can now start a linux environment using the WSL programm (ie. wsl.exe)

.. important:: Please note that all other commands in the manual should be executed in the WSL console 
  (and not in i.e. the windows powershell).