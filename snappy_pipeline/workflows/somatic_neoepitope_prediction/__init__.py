# -*- coding: utf-8 -*-
"""Implementation of the ``somatic_neoepitope_prediction`` step

The somatic_neoepitope_prediction step allows for the prediction of neoepitopes from somatic
(small) variant calling results and a transcript database such as ENSEMBL.  Further, the step
allows for the binding prediction to a given set of HLA alleles.

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

.. include:: DEFAULT_CONFIG_somatic_neoepitope_prediction.rst

"""

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

#: Default configuration for the somatic_gene_fusion_calling step
DEFAULT_CONFIG = r"""
step_config:
  somatic_neoepitope_prediction:
    path_somatic_variant_calling: REQUIRED  # REQUIRED
"""
