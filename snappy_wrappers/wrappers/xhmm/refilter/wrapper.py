# -*- coding: utf-8 -*-
# Merge the coverage tracks created in the "gatk_cov" steps.

from snakemake.shell import shell

shell(
    r"""
set -x

xhmm --matrix \
    -r {snakemake.input.original} \
    --excludeTargets {snakemake.input.filtered_targets_filter_center} \
    --excludeTargets {snakemake.input.filtered_targets_zscore_center} \
    --excludeSamples {snakemake.input.filtered_samples_filter_center} \
    --excludeSamples {snakemake.input.filtered_samples_zscore_center} \
    -o {snakemake.output}
"""
)
