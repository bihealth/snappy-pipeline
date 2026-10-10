# -*- coding: utf-8 -*-
"""Wrapper for running Defuse"""

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
    zcat {snakemake.input.reads_left} > $inputdir/reads_1.fastq.gz
fi
if [[ ! -f "$inputdir/reads_2.fastq.gz" ]]; then
    zcat {snakemake.input.reads_right} > $inputdir/reads_2.fastq.gz
fi

pushd $workdir

defuse_run.pl \
    -d {snakemake.input.dataset_directory} \
    -1 input/reads_1.fastq.gz \
    -2 input/reads_2.fastq.gz \
    -o output \
    -p 8

popd
"""
)
