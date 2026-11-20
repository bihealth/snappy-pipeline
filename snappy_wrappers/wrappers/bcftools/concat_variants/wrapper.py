# -*- coding: utf-8 -*-
"""Concatenate the VCFs of several tasks for one entity into one sorted, indexed VCF.

``bcftools concat`` needs the same samples in the same order in every input.
The calls of one pedigree or one library from different callers satisfy this.
"""

from typing import TYPE_CHECKING

from snappy_wrappers.snappy_wrapper import ShellWrapper

if TYPE_CHECKING:
    from snakemake.iocontainers import snakemake

ShellWrapper(snakemake).run(
    r"""
bcftools concat \
    --allow-overlaps \
    {extra_arguments} \
    {snakemake.input.vcf} \
    | bcftools sort \
        --temp-dir $TMPDIR \
        --output-type z \
        --output {snakemake.output.vcf}
tabix -f {snakemake.output.vcf}
"""
)
