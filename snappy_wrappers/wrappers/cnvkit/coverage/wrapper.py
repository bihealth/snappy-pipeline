# -*- coding: utf-8 -*-
"""Wrapper vor cnvkit.py coverage
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

cnvkit.py coverage {snakemake.input.bam} {snakemake.input.target}     -o {snakemake.output.target}
cnvkit.py coverage {snakemake.input.bam} {snakemake.input.antitarget} -o {snakemake.output.antitarget}

d=$(dirname "{snakemake.output.target}")
pushd $d
fn=$(basename "{snakemake.output.target}")
md5sum $fn > $fn.md5
fn=$(basename "{snakemake.output.antitarget}")
md5sum $fn > $fn.md5
popd
"""
)
