# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for external: Snakemake wrapper.py"""

from snappy_wrappers.snappy_wrapper import ShellWrapper

__author__ = "Oliver Stolpe <oliver.stolpe@bih-charite.de>"

args = getattr(snakemake.params, "args", {})

input = args["input"]
if not input:
    raise Exception("No bam found")

ShellWrapper(snakemake).run(
    r"""
set -x
mkdir -p $TMPDIR/tmp.d

# Link in bam files with the proper file name scheme
ln -sr {input} {snakemake.output.bam}

# Link in resultin BAM file or create index
if [[ -e {input}.bai ]]; then
    ln -sr {input}.bai {snakemake.output.bam_bai}
else
    samtools index {snakemake.output.bam}
fi

# QC Report ---------------------------------------------------------------------------------------

# gather statistics from BAM file
# TODO: use pipes for only reading once from disk?
samtools stats    {snakemake.output.bam} > {snakemake.output.report_bamstats_txt}
samtools flagstat {snakemake.output.bam} > {snakemake.output.report_flagstats_txt}
samtools idxstats {snakemake.output.bam} > {snakemake.output.report_idxstats_txt}


# Additional logging for transparency & reproducibility
# Logging: Save a copy this wrapper (with the pickle details in the header)
cp {__real_file__} $(dirname {snakemake.log.log})/wrapper.py

# Logging: Save a permanent copy of the environment file used
cp $(dirname {__file__})/environment.yaml $(dirname {snakemake.log.log})/environment_wrapper.yaml
"""
)
