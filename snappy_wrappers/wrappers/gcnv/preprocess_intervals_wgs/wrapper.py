# -*- coding: utf-8 -*-

from snakemake.shell import shell

shell(
    r"""
set -x

gatk PreprocessIntervals \
   --padding 0 \
   --interval-merging-rule OVERLAPPING_ONLY \
   --reference {snakemake.config[static_data_config][reference][path]} \
   --output {snakemake.output.interval_list}
"""
)
