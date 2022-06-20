# -*- coding: utf-8 -*-
"""Implementation of the ``ngs_sanity_checking`` step

Perform sanity checking from mapped reads for germline sample sheets, optionally taking the result
of ``hla_typing`` into consideration.

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

.. include:: DEFAULT_CONFIG_ngs_sanity_checking.rst

"""

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

#: Default configuration for the tcell_crg_report step
DEFAULT_CONFIG = r"""
# Default configuration ngs_sanity_checking
step_config:
  ngs_sanity_checking:
    path_ngs_mapping: ../path_ngs_mapping  # REQUIRED
    path_hla_typing: ../path_hla_typing    # OPTIONAl
    check_hla: true
"""
