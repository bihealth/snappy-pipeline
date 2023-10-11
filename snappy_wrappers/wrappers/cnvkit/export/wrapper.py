# -*- coding: utf-8 -*-
"""Wrapper for cnvkit.py export
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

set -x

export TMPDIR=$(mktemp -d)

cnvkit.py export seg {snakemake.input} -o {snakemake.output.seg}

cnvkit.py export vcf {snakemake.input} -o $TMPDIR/out.vcf

bgzip -c $TMPDIR/out.vcf > {snakemake.output.vcf}
tabix -f {snakemake.output.vcf}

rm -rf $TMPDIR

d=$(dirname "{snakemake.output.seg}")
pushd $d
fn=$(basename "{snakemake.output.seg}")
md5sum $fn > $fn.md5
fn=$(basename "{snakemake.output.vcf}")
md5sum $fn > $fn.md5
md5sum $fn.tbi > $fn.tbi.md5
popd
"""
)

# Compute MD5 sums of logs.
shell(
    r"""
md5sum {snakemake.log.log} >{snakemake.log.log_md5}
"""
)
