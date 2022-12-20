# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for target coverage report: Snakemake wrapper.py
"""

from snakemake import shell

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

shell.executable("/bin/bash")

shell(
    r"""
set -x

# Write out information about conda and save a copy of the wrapper with picked variables
# as well as the environment.yaml file.
conda list >{snakemake.log.conda_list}
conda info >{snakemake.log.conda_info}
md5sum {snakemake.log.conda_list} >{snakemake.log.conda_list_md5}
md5sum {snakemake.log.conda_info} >{snakemake.log.conda_info_md5}
cp {__real_file__} {snakemake.log.wrapper}
md5sum {snakemake.log.wrapper} >{snakemake.log.wrapper_md5}
cp $(dirname {__file__})/environment.yaml {snakemake.log.env_yaml}
md5sum {snakemake.log.env_yaml} >{snakemake.log.env_yaml_md5}

# Also pipe stderr to log file
if [[ -n "{snakemake.log.log}" ]]; then
    if [[ "$(set +e; tty; set -e)" != "" ]]; then
        rm -f "{snakemake.log.log}" && mkdir -p $(dirname {snakemake.log.log})
        exec 2> >(tee -a "{snakemake.log.log}" >&2)
    else
        rm -f "{snakemake.log.log}" && mkdir -p $(dirname {snakemake.log.log})
        echo "No tty, logging disabled" >"{snakemake.log.log}"
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

# Compute MD5 sums of logs.
shell(
    r"""
sleep 1s  # try to wait for log file flush
md5sum {snakemake.log.log} >{snakemake.log.log_md5}
"""
)
