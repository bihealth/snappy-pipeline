# -*- coding: utf-8 -*-
"""Wrapper for cnvkit.py coverage"""

from snakemake.shell import shell

__author__ = "Manuel Holtgrewe"
__email__ = "manuel.holtgrewe@bih-charite.de"

args = getattr(snakemake.params, "args", {})

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

# Function definitions ---------------------------------------------------------

coverage()
{{
    cnvkit.py coverage \
        --fasta {snakemake.input.reference} \
        --min-mapq {args[min_mapq]} \
        --processes {snakemake.threads} \
        {snakemake.input.bam} \
        --output $2 $1
}}

md5()
{{
    set -x

    fn=$1
    f=$(basename $fn)
    d=$(dirname $fn)
    pushd $d
    md5sum $f > $f.md5
    popd
}}

# -----------------------------------------------------------------------------

coverage {snakemake.input.target} {snakemake.output.target}
md5 {snakemake.output.target}

coverage {snakemake.input.antitarget} {snakemake.output.antitarget}
md5 {snakemake.output.antitarget}
"""
)

# Compute MD5 sums of logs.
shell(
    r"""
md5sum {snakemake.log.log} >{snakemake.log.log_md5}
"""
)
