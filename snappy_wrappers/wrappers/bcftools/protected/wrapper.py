# -*- coding: utf-8 -*-
"""Wrapper for running bcftools mpileup"""

from typing import TYPE_CHECKING

from snappy_wrappers.snappy_wrapper import ShellWrapper

if TYPE_CHECKING:
    from snakemake.iocontainers import snakemake

args = getattr(snakemake.params, "args", {})
filter_name = args["filter_name"]
bed = args["path_bed"]

ShellWrapper(snakemake).run(
    r"""
bcftools filter --soft-filter PROTECTED --mode + \
    --mask-file "{bed}" \
    -O z -o {snakemake.output.vcf} \
    {snakemake.input.vcf}
tabix {snakemake.output.vcf}
"""
)
