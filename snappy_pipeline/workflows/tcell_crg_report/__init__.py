# -*- coding: utf-8 -*-
"""Implementation of the ``tcell_crg_report`` step

This step collects all of the calls and information gathered for the BIH T cell CRG and generates
an Excel report for each patient.

.. note::

    Status: not implemented yet

==========
Step Input
==========

The BIH T cell CRG report generator uses the following as input:

- ``somatic_variant_annotation``
- ``somatic_epitope_prediction``
- ``somatic_variant_checking``
- ``somatic_ngs_sanity_checks``

===========
Step Output
===========

.. note:: TODO

=====================
Default Configuration
=====================

The default configuration is as follows.

.. include:: DEFAULT_CONFIG_tcell_crg_report.rst

=============================
Available Gene Fusion Callers
=============================

- ``cnvkit``

"""

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

#: Default configuration for the tcell_crg_report step
DEFAULT_CONFIG = r"""
# Default configuration tcell_crg_report
step_config:
  tcell_crg_report:
    path_somatic_variant_annotation: ../somatic_variant_annotation  # REQUIRED
"""
