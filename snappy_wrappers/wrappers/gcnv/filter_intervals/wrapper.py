# -*- coding: utf-8 -*-
from typing import TYPE_CHECKING

from snappy_wrappers.snappy_wrapper import ShellWrapper
from snappy_wrappers.tools.java import java_heap_mb

if TYPE_CHECKING:
    from snakemake.iocontainers import snakemake

args = snakemake.params.args
paths_cov = " ".join(snakemake.input.covs)

mem_jvm = java_heap_mb(snakemake)

ShellWrapper(snakemake).run(
    r"""
gatk --java-options "-Xmx{mem_jvm}m" \
    FilterIntervals \
    --intervals {snakemake.input.interval_list} \
    --annotated-intervals {snakemake.input.tsv} \
    --interval-merging-rule {args[interval_merging_rule]} \
    $(for tsv in {paths_cov}; do echo -I $tsv; done) \
    --minimum-gc-content {args[minimum_gc_content]} \
    --maximum-gc-content {args[maximum_gc_content]} \
    --minimum-mappability {args[minimum_mappability]} \
    --maximum-mappability {args[maximum_mappability]} \
    --output {snakemake.output.interval_list}
"""
)
