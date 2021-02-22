# -*- coding: utf-8 -*-
"""Implementation of the germline ``somatic_variant_checking`` step

The ``somatic_variant_checking`` step takes as the input the results of the
``somatic_variant_annotation`` step.  It then executes various tools computing statistics on the
result files and consistency checks with the pedigrees.

.. note::

    Status: not implemented yet

==========
Step Input
==========

The variant calling step uses Snakemake sub workflows for using the result of the
``somatic_variant_annotation`` step.

===========
Step Output
===========

.. note:: TODO

====================
Global Configuration
====================

.. note:: TODO

=====================
Default Configuration
=====================

The default configuration is as follows.

.. include:: DEFAULT_CONFIG_somatic_variant_checking.rst

=======
Reports
=======

Currently, no reports are generated.
"""

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>"

#: Available tools for checking variants
VARIANT_CHECKERS = ("bcftools_stats", "peddy")

#: Default configuration for the somatic_gene_fusion_calling step
DEFAULT_CONFIG = r"""
step_config:
  somatic_variant_checking:
    path_somatic_variant_calling: ../somatic_variant_calling  # REQUIRED
"""
