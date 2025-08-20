Troubleshooting errors
================
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

Reporting issues
----------------

If you encounter issues please:
 - make sure they are reproducible
 - indentify the step where the error occurs
 - rerun the workflow with an added "-vv" flag and record the output
   - i.e. ``stemcnv-check -vv > log.txt 2>&1``
 - Open a new issue on the `StemCNV-check GitHub repository <https://github.com/bihealth/StemCNV-check/>`_

 - You may need to create a GitHub account if you do not have one
 - make sure to attach the captured log output as well as the relevant log files produced by the snakemake rules (see 4.2)

