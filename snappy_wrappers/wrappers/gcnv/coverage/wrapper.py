# -*- coding: utf-8 -*-

import os
from snakemake.shell import shell

shell(
    r"""
set -x

gatk CollectReadCounts \
    --interval-merging-rule OVERLAPPING_ONLY \
    -R {snakemake.config[static_data_config][reference][path]} \
    -L {snakemake.input.interval_list} \
    -I {snakemake.input.bam} \
    --format TSV \
    -O {snakemake.output.tsv}
"""
)
