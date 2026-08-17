# -*- coding: utf-8 -*-
"""Wrapper for running bcftools mpileup"""

from typing import TYPE_CHECKING

from snappy_wrappers.snappy_wrapper import ShellWrapper

if TYPE_CHECKING:
    from snakemake.iocontainers import snakemake


args = getattr(snakemake.params, "args", {})
reference_path = args["reference_path"]
max_depth = args["max_depth"]

# FIXME: "locii" only ever gets set as the input, never as a parameter in args
if intervals := args["intervals"]:
    locii = "-r " + intervals
elif "locii" in snakemake.input.keys():
    locii = "-R " + snakemake.input.locii
elif locii_arg := args.get("locii"):
    locii = "-R " + locii_arg
else:
    locii = ""

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
