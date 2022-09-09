# -*- coding: utf-8 -*-
"""Wrapper for actually running ASCAT"""

from snakemake import shell

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

shell.executable("/bin/bash")

ARCAS_HLA_THREADS = 8
ARCAS_HLA_PAIRED_IS_PAIRED = True
_ARCAS_HLA_PAIRED = "--paired" if ARCAS_HLA_PAIRED_IS_PAIRED else ""

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

input=$(readlink -f {snakemake.input.bam})
mkdir -p work/star.arcashla.{snakemake.wildcards.library_name}/tmp/{{extracted,genotyped}}
mkdir -p $(dirname {snakemake.output.txt})
pushd work/star.arcashla.{snakemake.wildcards.library_name}

arcasHLA extract \
    {snakemake.input.bam} \
    -o tmp/extracted \
    {_ARCAS_HLA_PAIRED} \
    -t {ARCAS_HLA_THREADS} \
    -v

arcasHLA genotype \
    tmp/extracted/*.fq.gz \
    -o tmp/genotyped \
    -t {ARCAS_HLA_THREADS} \
    -v

popd

cp work/star.arcashla.{snakemake.wildcards.library_name}/tmp/genotyped/star.genotype.json {snakemake.output.txt}
pushd $(dirname {snakemake.output.txt})
md5sum $(basename {snakemake.output.txt}) >$(basename {snakemake.output.txt}).md5
"""
)
