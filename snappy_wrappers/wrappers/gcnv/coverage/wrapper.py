# -*- coding: utf-8 -*-

from snakemake.shell import shell

args = getattr(snakemake.params, "args", {})

shell(
    r"""
set -x

gatk CollectReadCounts \
    --interval-merging-rule OVERLAPPING_ONLY \
    -R {args[reference]} \
    -L {snakemake.input.interval_list} \
    -I {snakemake.input.bam} \
    --format TSV \
    -O {snakemake.output.tsv}
"""
)
