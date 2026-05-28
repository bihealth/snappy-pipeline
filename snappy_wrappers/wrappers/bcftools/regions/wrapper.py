# -*- coding: utf-8 -*-
"""Wrapper for running bcftools filter over regions defined by a bed file"""

from typing import TYPE_CHECKING

from snappy_wrappers.snappy_wrapper import ShellWrapper

if TYPE_CHECKING:
    from snakemake.iocontainers import snakemake

args = getattr(snakemake.params, "args", {})
filter_name = args["filter_name"]
bed = f"^{args['include']}" if "include" in args else args["exclude"]

ShellWrapper(snakemake).run(
    r"""
bcftools filter --soft-filter {filter_name} --mode + \
    --mask-file "{bed}" \
    -O z -o {snakemake.output.vcf} \
    {snakemake.input.vcf}
tabix {snakemake.output.vcf}
"""
)
