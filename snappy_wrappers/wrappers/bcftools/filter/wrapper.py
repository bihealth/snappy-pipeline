# -*- coding: utf-8 -*-
"""Wrapper for running bcftools filter"""

from typing import TYPE_CHECKING

from snappy_wrappers.snappy_wrapper import ShellWrapper

if TYPE_CHECKING:
    from snakemake.iocontainers import snakemake

args = getattr(snakemake.params, "args", {})
filter_name = args.get("filter_name", "bcftools")
mode = args.get("mode", "tag")
expression = (
    '--include "{}"'.format(args["include"])
    if args.get("include", None)
    else '--exclude "{}"'.format(args["exclude"])
)

if mode == "tag":
    cmd = r"""
bcftools filter --soft-filter {filter_name} --mode + \
    {expression} \
    -O z -o {snakemake.output.vcf} \
    {snakemake.input.vcf}
tabix {snakemake.output.vcf}
"""
else:
    cmd = r"""
bcftools filter \
    {expression} \
    -O z -o {snakemake.output.vcf} \
    {snakemake.input.vcf}
tabix {snakemake.output.vcf}
"""

ShellWrapper(snakemake).run(cmd)
