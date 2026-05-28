# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for Samtools - BAM QC report: Snakemake wrapper.py"""

from snakemake import shell

from snappy_wrappers.snappy_wrapper import ShellWrapper

shell.executable("/bin/bash")

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


# Build MD5 files for the reports
md5sum {snakemake.output.bamstats}  > {snakemake.output.bamstats_md5}
md5sum {snakemake.output.flagstats} > {snakemake.output.flagstats_md5}
md5sum {snakemake.output.idxstats}  > {snakemake.output.idxstats_md5}
"""
)

