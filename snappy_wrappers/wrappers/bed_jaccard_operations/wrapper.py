# -*- coding: utf-8 -*-
"""Wrapper for BEDTools intersect + snappy-bed_filter_jaccard: Snakemake wrapper.py
"""

from snakemake import shell

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

shell.executable("/bin/bash")

# Handle /dev/stdout output file
if str(snakemake.output) == ".wrapper.done":
    output_arg = ""
    output_bed = ""
else:
    output_arg = ">" + snakemake.output.bed
    output_bed = snakemake.output.bed

shell(
    r"""
set -x

# Also pipe stderr to log file
if [[ -n "{snakemake.log}" ]]; then
    if [[ "$(set +e; tty; set -e)" != "" ]]; then
        rm -f "{snakemake.log}" && mkdir -p $(dirname {snakemake.log})
        exec 2> >(tee -a "{snakemake.log}" >&2)
    else
        rm -f "{snakemake.log}" && mkdir -p $(dirname {snakemake.log})
        echo "No tty, logging disabled" >"{snakemake.log}"
    fi
fi

# Hack: get back bin directory of base/root environment.
export PATH=$PATH:$(dirname $(dirname $(which conda)))/bin

bedtools intersect \
    -wao \
    -a <(grep -v '^#' {snakemake.input.first} || true) \
    -b <({{ grep -v '^#' {snakemake.input.second} || true; }} | cut -f 1-3) \
| snappy-bed_filter_jaccard \
    --operation {snakemake.params.args[operation]} \
    --num-cols-first {snakemake.params.args[num_cols_first]} \
    --num-cols-second 3 \
    --threshold {snakemake.params.args[threshold]} \
| cut -f 1-{snakemake.params.args[num_cols_first]} \
{output_arg}

if [[ -n "{output_arg}" ]]; then
    pushd $(dirname {output_bed}) &&
        md5sum $(basename {output_bed}) >$(basename {output_bed}).md5
fi
"""
)
