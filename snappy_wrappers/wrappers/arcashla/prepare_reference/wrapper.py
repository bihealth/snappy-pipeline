# -*- coding: utf-8 -*-
"""Wrapper for preparing ARCAS-HLA reference"""

from typing import TYPE_CHECKING

from snappy_wrappers.snappy_wrapper import ShellWrapper

if TYPE_CHECKING:
    from snakemake.iocontainers import snakemake

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

ARCAS_HLA_REFERENCE_VERSION = "3.24.0"

ShellWrapper(snakemake).run(
    r"""
arcasHLA reference --version {ARCAS_HLA_REFERENCE_VERSION} --verbose

"""
)
