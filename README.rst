Coulter counter fitting scripts
===============================

The script attempts to fit both a mixture of normal, and a mixture of log normal
distributions for comparative purposes. It also displays plots adjusting for
the non-uniform bin width used by the Coulter counter.

Requirements
------------

To run the scripts, you'll need the following R packages:

* mixdist
* gridExtra
* tidyverse

Preparing the data
------------------

The script processes the CSV output files from a Coulter counter.  It will look
for any files with a ".CSV" extension in a subdirectory ``data`` relative to
the directory where the script is executed.

Running the script
------------------

The main script is ``fit_mixtures_to_coulter_data.R``. You should be able to
run it from the command line with:

.. code-block:: bash

    Rscript fit_mixtures_to_coulter_data.R

or interactively through R or Rstudio.

The script produces:

1. A file ``all_fittings.csv`` containing parameters of fitted log normal
   (prefixed with ``lnorm``) and normal (prefixed with ``norm``) distributions.
2. Per-input file plots showing volume fractions adjusted for bin width, with
   a ".pdf" extension.

