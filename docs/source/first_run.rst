Running the analysis for the first time
============

Before the first analysis sample table and config file need to be set up. Unless otherwise specified, stemcnv-check defaults to look for a "sample_table.tsv" (or .xlsx) and "config.yaml" file.

It is recommended to start by **creating a separate folder** for your project. This folder should include raw data folder, config.yaml and sample table files.

Once the config file is properly set up and the necessary static files are generated, you can run the StemCNV-check
workflow with simple command:
``stemcnv-check``

*Reminder for WSL users:* as before you may need to limit the memory usage of the workflow
and use this command or a variantion instead: ``stemcnv-check run -- --resources mem_mb=6000``
