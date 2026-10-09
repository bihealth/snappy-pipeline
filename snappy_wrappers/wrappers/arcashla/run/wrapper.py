# -*- coding: utf-8 -*-
"""Wrapper for running ARCAS-HLA"""

from typing import TYPE_CHECKING

from snappy_wrappers.snappy_wrapper import ShellWrapper

if TYPE_CHECKING:
    from snakemake.iocontainers import snakemake

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

ARCAS_HLA_THREADS = 8
ARCAS_HLA_PAIRED_IS_PAIRED = True
_ARCAS_HLA_PAIRED = "--paired" if ARCAS_HLA_PAIRED_IS_PAIRED else ""

ShellWrapper(snakemake).run(
    r"""
input=$(readlink -f {snakemake.input.bam})
mkdir -p work/{snakemake.wildcards.library_name}/tmp/{{extracted,genotyped}}
mkdir -p $(dirname {snakemake.output.txt})
pushd work/{snakemake.wildcards.library_name}/tmp

arcasHLA extract \
    {snakemake.input.bam} \
    -o extracted \
    {_ARCAS_HLA_PAIRED} \
    -t {ARCAS_HLA_THREADS} \
    -v

arcasHLA genotype \
    extracted/*.fq.gz \
    -o genotyped \
    -t {ARCAS_HLA_THREADS} \
    -v

popd

cp work/{snakemake.wildcards.library_name}/tmp/genotyped/*.genotype.json {snakemake.output.txt}
"""
)
