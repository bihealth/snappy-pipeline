# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for bcftools query: Snakemake wrapper.py"""

from typing import TYPE_CHECKING

from snappy_wrappers.snappy_wrapper import ShellWrapper

if TYPE_CHECKING:
    from snakemake.iocontainers import snakemake

__author__ = "Clemens Messerschmidt"

args = getattr(snakemake.params, "args", {})

ShellWrapper(snakemake).run(
    r"""
bcftools query \
    -f '{args[tumor_library]}\t%CHROM\t%POS\t%REF\t%ALT\n' \
    -s {args[tumor_library]} \
    {snakemake.input.vcf} \
> {snakemake.output.tsv}
"""
)
