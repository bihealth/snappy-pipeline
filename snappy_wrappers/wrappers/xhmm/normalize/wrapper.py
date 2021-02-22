# -*- coding: utf-8 -*-
# Perform normalization step after PCA.

from snakemake.shell import shell

shell(
    r"""
set -x

pca_dir=$(dirname {snakemake.input.pca})
pca_prefix=$(basename {snakemake.input.pca} .PC.txt)

xhmm \
    --normalize \
    -r {snakemake.input.centered} \
    --PCAfiles $pca_dir/$pca_prefix \
    --normalizeOutput {snakemake.output.normalized} \
    --PCnormalizeMethod PVE_mean \
    --PVE_mean_factor 0.7
"""
)
