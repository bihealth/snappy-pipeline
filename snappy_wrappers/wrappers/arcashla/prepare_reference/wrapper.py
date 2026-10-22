# -*- coding: utf-8 -*-
"""Wrapper for preparing ARCAS-HLA reference"""

from typing import TYPE_CHECKING

from snappy_wrappers.snappy_wrapper import ShellWrapper

if TYPE_CHECKING:
    from snakemake.iocontainers import snakemake

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

ShellWrapper(snakemake).run(
    r"""
arcasHLA reference --version {snakemake.params.reference_version} --verbose

"""
)
