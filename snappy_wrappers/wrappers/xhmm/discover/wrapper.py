# -*- coding: utf-8 -*-
# Perform variant discovery step for XHMM.

from snakemake.shell import shell

shell(
    r"""
set -x

out_dir=$(dirname {snakemake.output.xcnv})
out_prefix=$(basename {snakemake.output.xcnv} .xcnv)
tmp_dir=$out_dir/../tmp
mkdir -p $tmp_dir
echo -e "1e-08\t6\t70\t-3\t1\t0\t1\t3\t1" >$tmp_dir/params.txt

xhmm \
    --discover \
    -p $tmp_dir/params.txt \
    -r {snakemake.input.center_zscore} \
    -R {snakemake.input.refilter_original} \
    -c {snakemake.output.xcnv} \
    -a {snakemake.output.aux_xcnv} \
    -s $out_dir/$out_prefix
"""
)
