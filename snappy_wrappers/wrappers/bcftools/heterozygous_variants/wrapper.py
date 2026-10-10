# -*- coding: utf-8 -*-
"""Wrapper for finding heterozygous variants with bcftools"""

from typing import TYPE_CHECKING

from snappy_wrappers.snappy_wrapper import ShellWrapper

if TYPE_CHECKING:
    from snakemake.iocontainers import snakemake

args = getattr(snakemake.params, "args", {})

reference_path = snakemake.input.reference

# The tumor pileup is restricted to the heterozygous sites of the normal
locii = f"-R {snakemake.input.locii}" if "locii" in snakemake.input.keys() else ""

# Convert minimum B-allele fraction into ratio of alternative to reference alleles
min_ratio = args["min_baf"] / (1 - args["min_baf"])
max_ratio = 1 / min_ratio

min_depth = args["min_depth"]
max_depth = args["max_depth"]

ShellWrapper(snakemake).run(
    r"""
only_one_variant="N_ALT=2 & FORMAT/AD[:2]=0"
min_depth="FORMAT/AD[:0]>{min_depth} & FORMAT/AD[:1]>{min_depth}"
hetero="{min_ratio}*FORMAT/AD[:0]<=FORMAT/AD[:1] & FORMAT/AD[:1]<={max_ratio}*FORMAT/AD[:0]"

bcftools mpileup \
    {locii} \
    --max-depth {max_depth} \
    -f {reference_path} \
    -a "FORMAT/AD" \
    {snakemake.input.bam} \
    | bcftools filter \
    --include "$only_one_variant & $min_depth & $hetero" \
    -O z -o {snakemake.output.vcf}
tabix {snakemake.output.vcf}
"""
)
