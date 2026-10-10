# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for Samtools - BAM QC report: Snakemake wrapper.py"""

from snappy_wrappers.snappy_wrapper import ShellWrapper

ShellWrapper(snakemake).run(
    r"""
set -x

mkdir -p $TMPDIR/tmp.d

# QC Report
samtools stats    {snakemake.input.bam} > {snakemake.output.bamstats}
samtools flagstat {snakemake.input.bam} > {snakemake.output.flagstats}
samtools idxstats {snakemake.input.bam} > {snakemake.output.idxstats}
"""
)
