# -*- coding: utf-8 -*-
"""Wrapper for running bcftools filter"""

from typing import TYPE_CHECKING

from snappy_wrappers.snappy_wrapper import ShellWrapper

if TYPE_CHECKING:
    from snakemake.iocontainers import snakemake

args = getattr(snakemake.params, "args", {})
filter_name = args["filter_name"]
expression = (
    '--include "{}"'.format(args["include"])
    if args.get("include", None)
    else '--exclude "{}"'.format(args["exclude"])
)

ShellWrapper(snakemake).run(
    r"""
bcftools filter --soft-filter {filter_name} --mode + \
    {expression} \
    -O z -o {snakemake.output.vcf} \
    {snakemake.input.vcf}
tabix {snakemake.output.vcf}
"""
)
