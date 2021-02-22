# -*- coding: utf-8 -*-
# Merge the coverage tracks created in the "gatk_cov" steps.

from snakemake.shell import shell

shell(
    r"""
set -x

xhmm \
    --mergeGATKdepths \
    -o {snakemake.output} \
    $(for p in {snakemake.input}; do echo --GATKdepths $p; done)
"""
)
