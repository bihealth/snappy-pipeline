# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for FastQC: Snakemake wrapper.py"""

from typing import TYPE_CHECKING

from snappy_wrappers.snappy_wrapper import ShellWrapper

if TYPE_CHECKING:
    from snakemake.iocontainers import snakemake

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

args = getattr(snakemake.params, "args", {})

ShellWrapper(snakemake).run(
    r"""

outdir=$(dirname $(echo {snakemake.output}  | tr ' ' '\n' | tail -n 1))

_JAVA_OPTIONS="-Xms256m -Xmx512m -XX:CompressedClassSpaceSize=512m" \
fastqc \
    --noextract \
    -o $outdir \
    -t {args[num_threads]} \
    $(echo {snakemake.input} {args[more_reads]} | tr ' ' '\n' | grep 'fastq.gz$\|fastq$\|sam$\|bam$')

pushd $outdir
pwd
ls -lh
"""
)
