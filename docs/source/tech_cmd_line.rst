.. _tech-cmd-line:

Command line options
^^^^^^^^^^^^^^^^^^^^

.. admonition:: Under construction

    This page is still under construction and will contain more detailed information in the future

``stemcnv-check`` has 3 differemt subcommands:

  - ``setup-files``: Create example files for project setup
  - ``make-staticdata``: Create static data for SNP array
  - ``run``: Run the analysis workflow

Common command line options for all subcommands are:

.. option:: -h, --help

   Show this help message and exit

.. option:: --version

   Show version information and exit

.. option:: -v, --verbose

   More verbose output and additional logging, maximum verbosity at -vv

.. option:: -c, --config FILE

    Filename of config file. Default: config.yaml Note: if a global config exists in the cache path, it will also be used by default

.. option:: -s, --sample-table FILE

    Filename of sample table, can be tsv or xlsx format (1st sheet is read). Default: sample_table.tsv or sample_table.xlsx


Options specific to the ``setup-files`` subcommand are:


.. option:: --config-details {minimal,medium,advanced,complete}

    Level of detail for the config file. Default: minimal

.. option:: --sampletable-format {tsv,xlsx}

    Format of the sample table. Default: tsv

.. option:: --ovrrwrite

    Allow overwriting of existing files

Options specific to the ``make-staticdata`` *and* the ``run`` subcommands are:

.. option:: -d, --directory PATH

    Directory to run pipeline in. Default: current directory

.. option:: -n, --local-cores INT

    Number of cores for local submission. Default: 4

.. option:: --cache-path PATH

    Override auto-selection of the cache path to a specific directory. The default cache path is defined in the conifg file.

.. option:: --no-cache

    Do not use a cache directory for workflow created metadata. (cache includes: global array definition config, 
    conda envs, singularity images, and reference data). The default cache path is defined in the conifg file.

.. option:: --bind-points BIND_POINTS

    Additional bind points for apptainer containers, intended for expert users. Use i.e. '/path' to make it available 
    in apptainer, useful in case local directory contains symlinks that won't resolve in the container.


Options specific to only the ``run`` subcommand are:

.. option:: -t, --target {complete,report,collate-summary,summary-tables,collate-cnv-calls,combined-cnv-calls,PennCNV,CBS,SNP-data,gtc-data}

    Final target of the pipeline. Default: complete

.. option:: --collate-date [COLLATE_DATE]

    Add a date to the collate output files. Default without argument: today's date

.. option:: --snakemake-help

    Show snakemake help message, all snakemake options must be passed after '--'
