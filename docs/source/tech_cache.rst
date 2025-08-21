File caching
^^^^^^^^^^^^

StemCNV-check generates a cache folder (default: ``~/.cache/stemcnv-check``), where re-usable files that are not 
(necessarily) project specific are sac=ved for re-use. The location of the cache can be changed through the config file 
or with the ``--cache-path`` command line option.

The cache includes:
 - Genome reference files (fasta, gtf)
 - Array specific static data files (for easy re-use between projects)
 - conda environments generated and used by snakemake
 - apptainer images generated and used by snakemake 