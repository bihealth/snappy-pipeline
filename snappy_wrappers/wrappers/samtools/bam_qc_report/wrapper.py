# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for Samtools - BAM QC report: Snakemake wrapper.py"""

from snappy_wrappers.snappy_wrapper import ShellWrapper

args = getattr(snakemake.params, "args", {})

ShellWrapper(snakemake).run(
    r"""
set -x

mkdir -p $TMPDIR/tmp.d

# Validate input
if [[ "{args[bam_count]}" -eq  0 ]]; then
    echo "No BAM files provided!"
    exit 1
elif [[ "{args[bam_count]}" -gt  1 ]]; then
    echo "Multiple BAM files provided!"
    echo "{args[bam]}"
    exit 1
fi

# QC Report
samtools stats    {args[bam]} > {snakemake.output.bamstats}
samtools flagstat {args[bam]} > {snakemake.output.flagstats}
samtools idxstats {args[bam]} > {snakemake.output.idxstats}
"""
)

