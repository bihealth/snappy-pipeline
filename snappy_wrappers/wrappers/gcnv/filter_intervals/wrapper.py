# -*- coding: utf-8 -*-
from math import ceil
from typing import TYPE_CHECKING

from snappy_wrappers.snappy_wrapper import ShellWrapper

if TYPE_CHECKING:
    from snakemake.iocontainers import snakemake

paths_cov = " ".join(snakemake.input.covs)

# As an estimate, give 75% of total memory as usable to JVM
mem_jvm = ceil(snakemake.resources.mem_mb * 0.75)

ShellWrapper(snakemake).run(
    r"""

gatk --java-options "-Xmx{mem_jvm}m" \
    FilterIntervals \
    --intervals {snakemake.input.interval_list} \
    --annotated-intervals {snakemake.input.tsv} \
    --interval-merging-rule OVERLAPPING_ONLY \
    $(for tsv in {paths_cov}; do echo -I $tsv; done) \
    --minimum-gc-content 0.1 \
    --maximum-gc-content 0.9 \
    --minimum-mappability 0.9 \
    --maximum-mappability 1.0 \
    --output {snakemake.output.interval_list}
"""
)
