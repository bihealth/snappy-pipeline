# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for bcftools stats: Snakemake wrapper.py
"""

from snakemake import shell

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

shell(
    r"""
set -x

export TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT

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

bcftools stats \
    --split-by-ID \
    -s {snakemake.wildcards.donor_ngs_library} \
    {snakemake.input.vcf} \
> $(dirname {snakemake.output.txt})/$(basename {snakemake.output.txt} .txt).txt

pushd $(dirname {snakemake.output.txt}) && \
    md5sum $(basename {snakemake.output.txt} .txt).txt \
    > $(basename {snakemake.output.txt} .txt).txt.md5 && \
    popd
"""
)
