Troubleshooting errors


Common issues running StemCNV-check
------------

The StemCNV-check tool runs in two parts. 

First it parses the sample table and config files and prepares the actual analysis,
which is then started for the second part, which internally uses snakemake to manage the workflow.
Once the second part starts snakemake will first write "Building DAG of jobs..." and the continuous progress with individual 
steps to the terminal.

- **Errors before workflow execution by snakemake**

The most likely source of errors before the workflow starts are formatting issues or values outside allowed configuration
in the sample table or config file. In these cases StemCNV-check should tell you what problems exist, so that you can fix them.
Please read these error messages carefully and adjust the files accordingly.  
Any error messages including or originating from the yaml (parser) packages also indicate a problem with the config file.


- **Errors during workflow execution by snakemake**

If the workflow fails during execution, snakemake automatically reports which step failed. If this happens there is most 
likely an unforseen issue with the workflow. To help us fix such issues, please make sure to include the corresponding 
log files (indicating by snakemake in the step that failed) when reporting the issue. Since such errors are often specific 
to the data it may be necessary to provide the raw data and config files to reproduce the error.


-**Snakemake: Directory cannot be locked**

If you encounter an error message like "Directory cannot be locked" then you either attempted to start the workflow 
twice in the same directory (at the same time), or a previous run of the workflow was not properly finished. 
The latter can be caused by shutting down the running workflow, either from a computer shutdown or from closing the 
terminal window, in which the workflow was running. In this case the snakemake workflow manager will refuse to restart 
until the improperly removed previous lock is removed. This can be done using this command: 
``stemcnv-check run -- --unlock``

If you need to shut down yur computer but the workflow ist still running you can interrupt it by using ``Ctrl+C`` key 
combination in the respective terminal window. On Linux this means "keyboard interrupt" and will instruct the workflow 
to abort. The abortion process is not instantaneous and may take a few moments to stop all running processes.

-**WSL: Unclear crashes & memory overload**

If StemCNV-check crashes unexpectedly on WSL or the system becomes unresponsive, this can be caused by too much memory 
usage inside the virtual WSL system. Rules exiting with a SIGKILL signal are also an indicator of this.  
To prevent this it is advisable to specifically instruct StemCNV-check how much memory can be used (otherwise it will 
use all available, which may not be properly defined or detected inside WSL). Memory should be given in Mb, so for i.e.
6Gb of memory use:  
``stemcnv-check {--options} run -- --resources mem_mb=6000``


-**Downloads / Internet connection not working**

StemCNV-check requires an internet connection to download necessary files (both genome references and package environments
from conda). If you have a restricted internet connection, you may need to configure a proxy server (see below). 
An alternative solution to this is to connect to an open network (i.e. eduroam which is widely available) and run the 
`stemcnv-check make-staticdata` command and at least one (example) run of the `stemcnv-check [run]` workflow. 
After this, all necessary files should be available locally and the internet connection is no longer required.

If your institute has implemented a proxy server to facilitate unrestricted internet connections, follow the respective 
instructions.  

Windows proxy settings are generally not automatically applied to WSL. 
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


You can make these changes either manually by using a text editor like nano (`nano ~/.bashrc`) or by using the following commands:

.. code:: bash

    echo "export http_proxy=http://{proxy.your_institute.com}:{port}" >> ~/.bashrc
    echo "export https_proxy=http://{proxy.your_institute.com}:{port}" >> ~/.bashrc
    echo "export ftp_proxy=http://{proxy.your_institute.com}:{port}" >> ~/.bashrc
    echo "export PIP_PROXY=http://{proxy.your_institute.com}:{port}" >> ~/.bashrc
    echo "export no_proxy=$no_proxy,localhost,{intranet.domain,your_institute.com}" >> ~/.bashrc



-**Reporting issues**

If you encounter issues please:
 - make sure they are reproducible
 - indentify the step where the error occurs
 - rerun the workflow with an added "-vv" flag and record the output
   - i.e. ``stemcnv-check -vv > log.txt 2>&1``
 - Open a new issue on the `StemCNV-check GitHub repository<https://github.com/bihealth/StemCNV-check/issues>`_
   - You may need to create a GitHub account if you do not have one
 - make sure to attach the captured log output as well as the relevant log files produced by the snakemake rules (see 4.2)

