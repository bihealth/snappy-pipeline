# -*- coding: utf-8 -*-

from typing import TYPE_CHECKING

from snappy_wrappers.snappy_wrapper import ShellWrapper

if TYPE_CHECKING:
    from snakemake.iocontainers import snakemake

args = getattr(snakemake.params, "args", {})

ShellWrapper(snakemake).run(
    r"""

gatk CollectReadCounts \
    --interval-merging-rule OVERLAPPING_ONLY \
    -R {args[reference]} \
    -L {snakemake.input.interval_list} \
    -I {snakemake.input.bam} \
    --format TSV \
    -O {snakemake.output.tsv}
"""
)
