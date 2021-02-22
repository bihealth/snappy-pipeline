# -*- coding: utf-8 -*-
# Perform the PCA step of XHMM.

from snakemake.shell import shell

shell(
    r"""
set -x

out_dir=$(dirname {snakemake.output.pc})
out_prefix=$(basename {snakemake.output.pc} .PC.txt)

xhmm \
    --PCA \
    -r {snakemake.input} \
    --PCAfiles $out_dir/$out_prefix
"""
)
