# -*- coding: utf-8 -*-
"""Gather the VCF chunks of a scatter-gather run into one sorted, indexed VCF."""

from typing import TYPE_CHECKING

from snappy_wrappers.snappy_wrapper import ShellWrapper

if TYPE_CHECKING:
    from snakemake.iocontainers import snakemake

ShellWrapper(snakemake).run(
    r"""
bcftools concat \
    --allow-overlaps \
    --rm-dups none \
    {snakemake.input.vcf} \
| bcftools sort \
    --temp-dir $TMPDIR \
    --output-type z \
    --output {snakemake.output.vcf}
tabix -f {snakemake.output.vcf}
"""
)
