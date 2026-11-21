# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for GetPileupSummaries: Snakemake wrapper.py"""

from snappy_wrappers.snappy_wrapper import ShellWrapper

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

args = getattr(snakemake.params, "args", {})

reference = snakemake.input.reference
common_variants = snakemake.input.common_variants

if java_options := args.get("java_options", None):
    java_options = f"--java-options '{java_options}'"

ShellWrapper(snakemake).run(
    r"""
# export JAVA_HOME=$(dirname $(which gatk))/..
export LD_LIBRARY_PATH=$(dirname $(which bgzip))/../lib

gatk {java_options} GetPileupSummaries \
    --input {snakemake.input.bam} \
    --reference {reference} \
    --variant {common_variants} \
    --intervals {common_variants} \
    --output {snakemake.output.pileup} \
    {extra_arguments}
"""
)
