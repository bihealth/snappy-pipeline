# -*- coding: utf-8 -*-
from math import ceil

from snakemake.shell import shell

paths_cov = " ".join(snakemake.input.covs)

# Extract resource usage in megabytes
mem = str(snakemake.resources.memory).lower()
for i, suffix in enumerate("kmg"):
    if mem.endswith(suffix):
        mem_i = int(mem[:-1]) * (2 ** (10 * (i + 1)))
        break
else:
    mem_i = int(mem)

# As an estimate, give 75% of total memory as usable to JVM
mem_jvm = ceil(mem_i / 1024.0 / 1024.0 * 0.75)

shell(
    r"""
set -x



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
