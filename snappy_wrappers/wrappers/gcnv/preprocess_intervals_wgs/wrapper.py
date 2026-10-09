# -*- coding: utf-8 -*-

from typing import TYPE_CHECKING

from snappy_wrappers.snappy_wrapper import ShellWrapper

if TYPE_CHECKING:
    from snakemake.iocontainers import snakemake

ShellWrapper(snakemake).run(
    r"""

gatk PreprocessIntervals \
   --padding 0 \
   --interval-merging-rule OVERLAPPING_ONLY \
   --reference {snakemake.input.reference} \
   --output {snakemake.output.interval_list}
"""
)
