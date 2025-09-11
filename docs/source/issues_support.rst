.. _issues-support:

Reporting issues
^^^^^^^^^^^^^^^^

.. admonition:: Under construction

    This page is still under construction and will contain more detailed information in the future
    

If you encounter issues please:
 - Rerun the worfklow to make sure your issues are reproducible
 - Indentify the step where the error occurs (see :ref:`StemCNV-check error descriptions <issues-intro>`)
 - Rerun the workflow with an added "-vv" flag and record the output

  .. code-block:: bash
  
    stemcnv-check run -vv > log.txt 2>&1
 
 - Open a new issue on the `StemCNV-check GitHub repository <https://github.com/bihealth/StemCNV-check/>`_
    - You may need to create a GitHub account if you do not have one
 - Make sure to attach the captured `log.txt` output as well as the relevant log files produced by the snakemake rules

