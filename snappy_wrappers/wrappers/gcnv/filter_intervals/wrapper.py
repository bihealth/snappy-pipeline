# -*- coding: utf-8 -*-
from snakemake.shell import shell

paths_cov = " ".join(snakemake.input.covs)

shell(
    r"""
set -x

gatk FilterIntervals \
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
