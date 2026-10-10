# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for Hera: Snakemake wrapper.py"""

from typing import TYPE_CHECKING

from snappy_wrappers.snappy_wrapper import ShellWrapper

if TYPE_CHECKING:
    from snakemake.iocontainers import snakemake

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

ShellWrapper(snakemake).run(
    r"""
echo ${{JOB_ID:-unknown}} >$(dirname {snakemake.output.done})/sge_job_id


workdir=$(dirname {snakemake.output.done})
inputdir=$workdir/input

mkdir -p $inputdir

if [[ ! -f "$inputdir/reads_1.fastq.gz" ]]; then
    cat {snakemake.input.reads_left} > $inputdir/reads_1.fastq.gz
fi
if [[ ! -f "$inputdir/reads_2.fastq.gz" ]]; then
    cat {snakemake.input.reads_right} > $inputdir/reads_2.fastq.gz
fi

pushd $workdir

hera quant \
    -i {snakemake.input.index} \
    -f {snakemake.input.genome} \
    -t 8 \
    -o $PWD \
    input/reads_1.fastq.gz \
    input/reads_2.fastq.gz

popd
"""
)
