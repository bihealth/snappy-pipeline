# -*- coding: utf-8 -*-
"""Implementation of the ``somatic_variant_expression`` step

This step allows the combination of somatic variant calling results with their expression
from RNA-seq data.  This allows for (1) extending a somatic VCF file with columns for the
corresponding RNA-seq data giving depth of coverage and minor allele fraction in the tumor
RNA-eq and (2) for computing a p value for likelihood of observation by chance.

.. note::

    Status: not implemented yet

==========
Step Input
==========

.. note:: TODO

===========
Step Output
===========

.. note:: TODO

=====================
Default Configuration
=====================

The default configuration is as follows.

.. include:: DEFAULT_CONFIG_somatic_variant_expression.rst

"""

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>"

#: Default configuration for the somatic_variant_expression step
DEFAULT_CONFIG = r"""
step_config:
  somatic_variant_expression:
    path_ngs_mapping: ../ngs_mapping                          # REQUIRED
    path_somatic_variant_calling: ../somatic_variant_calling  # REQUIRED
"""
