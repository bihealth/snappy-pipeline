# -*- coding: utf-8 -*-

from typing import TYPE_CHECKING

from snappy_wrappers.snappy_wrapper import ShellWrapper

if TYPE_CHECKING:
    from snakemake.iocontainers import snakemake

args = snakemake.params.args

ShellWrapper(snakemake).run(
    r"""

gatk PreprocessIntervals \
   --padding {args[padding]} \
   --interval-merging-rule {args[interval_merging_rule]} \
   --reference {snakemake.input.reference} \
   --output {snakemake.output.interval_list}
"""
)
