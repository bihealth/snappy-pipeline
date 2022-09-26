# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for target coverage report: Snakemake wrapper.py
"""

from snakemake import shell

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

shell.executable("/bin/bash")


shell(
    r"""
set -x

# Hack: get back bin directory of base/root environment.
export PATH=$PATH:$(dirname $(dirname $(which conda)))/bin

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

export TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT

# Get sorted targets BED file.
zcat --force {snakemake.params.args[path_targets_bed]} \
| awk -F $'\t' 'BEGIN {{ OFS = FS; }} ($2 < $3) {{ print; }}' \
> $TMPDIR/targets.tmp.bed

bedtools sort \
    -i $TMPDIR/targets.tmp.bed \
    -faidx {snakemake.config[static_data_config][reference][path]}.genome \
| uniq \
> $TMPDIR/targets.bed

bedtools coverage \
    -a $TMPDIR/targets.bed \
    -b {snakemake.input.bam} \
    -g {snakemake.config[static_data_config][reference][path]}.genome \
    -hist \
    -sorted \
| python $(dirname {__file__})/../../../tools/bam_cov_stats.py \
    --bed-path $TMPDIR/targets.bed \
    --min-cov-warning {snakemake.params.args[min_cov_warning]} \
    --min-cov-ok {snakemake.params.args[min_cov_ok]} \
    --max-coverage {snakemake.params.args[max_coverage]} \
    $(if [[ "{snakemake.params.args[detailed_reporting]}" == "True" ]]; then \
        echo --report dec; \
    fi) \
> {snakemake.output.txt}

md5sum {snakemake.output.txt} > {snakemake.output.txt_md5}
"""
)
