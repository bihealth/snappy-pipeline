# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for bcftools stats: Snakemake wrapper.py"""

from snakemake import shell

__author__ = "Clemens Messerschmidt"

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

bcftools query \
    -f '{snakemake.wildcards.tumor_library}\t%CHROM\t%POS\t%REF\t%ALT\n' \
    -s {snakemake.wildcards.tumor_library} \
    {snakemake.input.vcf} \
> {snakemake.output.tsv}

md5sum {snakemake.output.tsv} > {snakemake.output.tsv}.md5
"""
)
