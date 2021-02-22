# -*- coding: utf-8 -*-
# Perform Z-score computation and centering.

from snakemake.shell import shell

shell(
    r"""
set -x

xhmm \
    --matrix \
    -r {snakemake.input} \
    --centerData \
    --centerType sample \
    --zScoreData \
    -o {snakemake.output.zscore_center} \
    --outputExcludedTargets {snakemake.output.filtered_targets} \
    --outputExcludedSamples {snakemake.output.filtered_samples} \
    --maxSdTargetRD 30
"""
)
