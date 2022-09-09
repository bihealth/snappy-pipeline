# -*- coding: utf-8 -*-
"""Wrapper vor cnvkit.py report
"""

from snakemake.shell import shell

__author__ = "Manuel Holtgrewe"
__email__ = "manuel.holtgrewe@bih-charite.de"

shell(
    r"""
# Also pipe everything to log file
if [[ -n "{snakemake.log.log}" ]]; then
    if [[ "$(set +e; tty; set -e)" != "" ]]; then
        rm -f "{snakemake.log.log}" && mkdir -p $(dirname {snakemake.log.log})
        exec &> >(tee -a "{snakemake.log.log}" >&2)
    else
        rm -f "{snakemake.log.log}" && mkdir -p $(dirname {snakemake.log.log})
        echo "No tty, logging disabled" >"{snakemake.log.log}"
    fi
fi

# Write out information about conda installation.
conda list >{snakemake.log.conda_list}
conda info >{snakemake.log.conda_info}
md5sum {snakemake.log.conda_list} >{snakemake.log.conda_list_md5}
md5sum {snakemake.log.conda_info} >{snakemake.log.conda_info_md5}

# -----------------------------------------------------------------------------

cnvkit.py breaks {snakemake.input.cnr} {snakemake.input.cns} \
> {snakemake.output.breaks}

cnvkit.py gainloss {snakemake.input.cnr} -s {snakemake.input.cns} \
> {snakemake.output.gainloss}

cnvkit.py gender {snakemake.input.cnr} \
> {snakemake.output.gender}

cnvkit.py metrics {snakemake.input.cnr} -s {snakemake.input.cns} \
> {snakemake.output.metrics}

cnvkit.py segmetrics {snakemake.input.cnr} \
    --mean \
    --median \
    --mode \
    --stdev \
    --sem \
    --mad \
    --mse \
    --iqr \
    --bivar \
    --ci \
    --pi \
    -s {snakemake.input.cns} \
    -o {snakemake.output.segmetrics}

d=$(dirname "{snakemake.output.breaks}")
pushd $d
fn=$(basename "{snakemake.output.breaks}")
md5sum $fn > $fn.md5
fn=$(basename "{snakemake.output.gainloss}")
md5sum $fn > $fn.md5
fn=$(basename "{snakemake.output.gender}")
md5sum $fn > $fn.md5
fn=$(basename "{snakemake.output.metrics}")
md5sum $fn > $fn.md5
fn=$(basename "{snakemake.output.segmetrics}")
md5sum $fn > $fn.md5
popd
"""
)
