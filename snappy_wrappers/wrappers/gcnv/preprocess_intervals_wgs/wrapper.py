# -*- coding: utf-8 -*-

from typing import TYPE_CHECKING

from snappy_wrappers.snappy_wrapper import ShellWrapper

if TYPE_CHECKING:
    from snakemake.iocontainers import snakemake

args = getattr(snakemake.params, "args", {})

ShellWrapper(snakemake).run(
    r"""

gatk PreprocessIntervals \
   --padding 0 \
   --interval-merging-rule OVERLAPPING_ONLY \
   --reference {args[reference]} \
   --output {snakemake.output.interval_list}
"""
)
