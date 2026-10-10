# -*- coding: utf-8 -*-
"""Wrapper for running bcftools mpileup"""

from typing import TYPE_CHECKING

from snappy_wrappers.snappy_wrapper import ShellWrapper

if TYPE_CHECKING:
    from snakemake.iocontainers import snakemake


args = getattr(snakemake.params, "args", {})
reference_path = snakemake.input.reference
max_depth = args["max_depth"]

# The tumor pileup is restricted to the heterozygous sites of the normal
locii = f"-R {snakemake.input.locii}" if "locii" in snakemake.input.keys() else ""

ShellWrapper(snakemake).run(
    r"""
bcftools mpileup \
    {locii} \
    --max-depth {max_depth} \
    -f {reference_path} \
    -a "FORMAT/AD" \
    -O z -o {snakemake.output.vcf} \
    {snakemake.input.bam}
tabix {snakemake.output.vcf}
"""
)
