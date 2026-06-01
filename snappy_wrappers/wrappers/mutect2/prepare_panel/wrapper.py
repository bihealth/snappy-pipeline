# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for MuTect 2: Snakemake wrapper.py"""

from snappy_wrappers.snappy_wrapper import ShellWrapper

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

args = getattr(snakemake.params, "args", {})

intervals = getattr(snakemake.input, "intervals", "")
if intervals:
    intervals = f"--intervals {intervals}"

if java_options := args.get("java_options", ""):
    java_options = f"--java-options '{java_options}'"

extra_arguments = " ".join(args.get("extra_arguments", []))

ShellWrapper(snakemake).run(
    r"""
set -x

export JAVA_HOME=$(dirname $(which gatk))/..
export LD_LIBRARY_PATH=$(dirname $(which bgzip))/../lib

vcf=$(basename --suffix=.gz {snakemake.output.vcf})

gatk {java_options} Mutect2 \
    --tmp-dir $TMPDIR \
    --reference {snakemake.input.reference} \
    --input {snakemake.input.normal_bam} \
    {intervals} \
    --max-mnp-distance 0 \
    --output $TMPDIR/$vcf \
    {extra_arguments}

bgzip $TMPDIR/$vcf
tabix $TMPDIR/$vcf.gz

mv $TMPDIR/$vcf.gz $TMPDIR/$vcf.gz.tbi $(dirname {snakemake.output.vcf})
"""
)
