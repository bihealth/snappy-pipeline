# -*- coding: utf-8 -*-
# Use GATK 3 for filtering and centering the merged coverages.

from snakemake.shell import shell

shell(
    r"""
set -x

xhmm \
    --matrix \
    -r {snakemake.input.merge_cov} \
    --centerData \
    --centerType target \
    -o {snakemake.output.centered} \
    --outputExcludedTargets {snakemake.output.filtered_targets} \
    --outputExcludedSamples {snakemake.output.filtered_samples} \
    --excludeTargets {snakemake.input.extreme_gc} \
    --minTargetSize 10 \
    --maxTargetSize 10000 \
    --minMeanTargetRD 10 \
    --maxMeanTargetRD 500 \
    --minMeanSampleRD 25 \
    --maxMeanSampleRD 200 \
    --maxSdSampleRD 150
"""
)
