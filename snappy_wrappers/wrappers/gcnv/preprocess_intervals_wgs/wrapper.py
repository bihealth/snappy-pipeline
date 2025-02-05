# -*- coding: utf-8 -*-

from snakemake.shell import shell

args = getattr(snakemake.params, "args", {})

shell(
    r"""
set -x

gatk PreprocessIntervals \
   --padding 0 \
   --interval-merging-rule OVERLAPPING_ONLY \
   --reference {args[reference]} \
   --output {snakemake.output.interval_list}
"""
)
