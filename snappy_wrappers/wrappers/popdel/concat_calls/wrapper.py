# -*- coding: utf-8 -*-
"""Wrapper for concatenating the chromosome-wise "popdel call" output files."""

from typing import TYPE_CHECKING

from snappy_wrappers.snappy_wrapper import ShellWrapper

if TYPE_CHECKING:
    from snakemake.iocontainers import snakemake

__author__ = "Manuel Holtgrewe"
__email__ = "manuel.holtgrewe@bih-charite.de"

ShellWrapper(snakemake).run(
    r"""
bcftools concat {snakemake.input.vcf} \
| bcftools sort -O z -o {snakemake.output.vcf}

tabix -f {snakemake.output.vcf}
"""
)
