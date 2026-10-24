# -*- coding: utf-8 -*-
"""Wrapper for running STAR-Fusion"""

from typing import TYPE_CHECKING

from snappy_wrappers.snappy_wrapper import ShellWrapper

if TYPE_CHECKING:
    from snakemake.iocontainers import snakemake

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

# The reads are inputs. If another task wrote them, the input is that task's out directory and
# params list the files.
args = getattr(snakemake.params, "args", {})
reads = args.get("input", {})
reads_left = snakemake.input.get("reads_left") or reads["reads_left"]
reads_right = snakemake.input.get("reads_right") or reads.get("reads_right", "")

ShellWrapper(snakemake).run(
    r"""
workdir={snakemake.output.out}
mkdir -p $workdir
echo ${{JOB_ID:-unknown}} > $workdir/sge_job_id

inputdir=$workdir/input

mkdir -p $inputdir

if [[ ! -f "$inputdir/reads_1.fastq.gz" ]]; then
    cat {reads_left} > $inputdir/reads_1.fastq.gz
fi
if [[ ! -f "$inputdir/reads_2.fastq.gz" ]]; then
    cat {reads_right} > $inputdir/reads_2.fastq.gz
fi

pushd $workdir

mkdir -p output

STAR-Fusion \
    --genome_lib_dir {snakemake.input.ctat_resource_lib} \
    --left_fq input/reads_1.fastq.gz \
    --right_fq input/reads_2.fastq.gz \
    --output_dir output \
    {extra_arguments}

popd
"""
)
