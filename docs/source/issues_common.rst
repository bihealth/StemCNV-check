.. _issues-common:

Common issues and solutions
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Snakemake: Directory cannot be locked
=====================================

If you encounter an error message like "Directory cannot be locked" then you either attempted to start the workflow 
twice in the same directory (at the same time), or a previous run of the workflow was not properly finished. 
The latter can be caused by shutting down the running workflow, either from a computer shutdown or from closing the 
terminal window, in which the workflow was running. In this case the snakemake workflow manager will refuse to restart 
until the improperly removed previous lock is removed. This can be done using this command: 

.. code-block:: bash

    stemcnv-check run -- --unlock

If you need to shut down your computer but the workflow ist still running you can interrupt it by using ``Ctrl+C`` key 
combination in the respective terminal window. On Linux this means "keyboard interrupt" and will instruct the workflow 
to abort. The abortion process is not instantaneous and may take a few moments to stop all running processes.

Unclear crashes & memory overload in a VM or WSL
================================================

To run StemCNV-check on windows you may be using a virtual machine (or WSL) to provide a linux environment. 
In this case, StemCNV-check running inside the (virtual) linux may not always be aware of how much memory is actually 
available on the computer, which can then lead to memory over-usage and crashes.  
On WSL, StemCNV-check should automatically detect the available memory and try to restrict itself to that, however this 
is not an exact process and some over-usage may still occur.

If you observe frequent and unexpected crashes when running StemCNV-check on a VM or WSL, especially when also the host 
system becomes unresponsive, this may be caused by too much memory usage inside the virtual WSL system. Rules exiting 
with a SIGKILL signal are also an indicator of this.  

To prevent this you can instruct StemCNV-check how much memory can be used (the default is to sue all available memory, 
which should usually work well on physical linux or WSL). This can be done with the ``--memory-mb`` option, which takes 
a number (in Mb) for the total amount of memory to use. Memory usage in StemNV-check (and the underlying snakemake) is 
primarily estimated, so you may want to allow for some buffer (i.e. if your VM can access 8GB of memory, only give 7-7.5 
to StemCNV-check.). An example command would then be:

.. code-block:: bash

    stemcnv-check run --memory-mb=7500


Downloads / Internet connection not working
===========================================

StemCNV-check requires an internet connection to download necessary files (both genome references and package environments
from conda). If you have a restricted internet connection, you may need to configure a proxy server (see below). 
An alternative solution to this is to connect to an open network (i.e. eduroam which is widely available) and run the 
``stemcnv-check make-staticdata`` command and at least one (example) run of the ``stemcnv-check run`` workflow. 
After this, all necessary files should be available locally and the internet connection is no longer required.

If your institute has implemented a proxy server to facilitate unrestricted internet connections, follow the respective 
instructions.  

(Automatic) Windows proxy settings are not always automatically applied to WSL. 
Therefore, you may need to manually configure WSL to use the proxy. To do this, first identify the proxy server 
(i.e. "proxy.your_institute.com") and port (usually "8080") that should be used in your network. Then add the 
following lines to your .bashrc file in your WSL home directory (i.e. /home/username/.bashrc)

.. code:: bash

    export http_proxy=http://{proxy.your_institute.com}:{port}
    export https_proxy=http://{proxy.your_institute.com}:{port}
    export ftp_proxy=http://{proxy.your_institute.com}:{port}
    export PIP_PROXY=http://{proxy.your_institute.com}:{port}


You may also want to exclude certain (intranet) addresses that are reachable through your normal connection, 
but possibly not through the proxy that allows internet access (i.e. "your_institute.com"). Multiple addresses can be 
entered if they are separated by commas. This can be done with the following line:

.. code:: bash

    export no_proxy=$no_proxy,localhost,{intranet.domain,your_institute.com}


You can make these changes either manually by using a text editor like nano (``nano ~/.bashrc``) or by using the following commands:

.. code:: bash

    echo "export http_proxy=http://{proxy.your_institute.com}:{port}" >> ~/.bashrc
    echo "export https_proxy=http://{proxy.your_institute.com}:{port}" >> ~/.bashrc
    echo "export ftp_proxy=http://{proxy.your_institute.com}:{port}" >> ~/.bashrc
    echo "export PIP_PROXY=http://{proxy.your_institute.com}:{port}" >> ~/.bashrc
    echo "export no_proxy=$no_proxy,localhost,{intranet.domain,your_institute.com}" >> ~/.bashrc
