# -*- coding: utf-8 -*-

from typing import TYPE_CHECKING

from snappy_wrappers.snappy_wrapper import ShellWrapper

if TYPE_CHECKING:
    from snakemake.iocontainers import snakemake

ShellWrapper(snakemake).run(
    r"""

gatk CollectReadCounts \
    --interval-merging-rule OVERLAPPING_ONLY \
    -R {snakemake.input.reference} \
    -L {snakemake.input.interval_list} \
    -I {snakemake.input.bam} \
    --format TSV \
    -O {snakemake.output.tsv}
"""
)
