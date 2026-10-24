# -*- coding: utf-8 -*-
"""Wrapper to select the tumor (and normal) libraries of a vcf, rename them to samples & split multiallelics"""

from snappy_wrappers.snappy_wrapper import ShellWrapper

__author__ = "Eric Blanc"
__email__ = "eric.blanc@bih-charite.de"

args = getattr(snakemake.params, "args", {})

# Renaming libraries to samples, and ensure normal before tumor
# (bcftools view ensures order of libraries, bcftools reheader renames)
samples = args["tumor_sample"]
libraries = args["tumor_library"]
if "normal_sample" in args and "normal_library" in args:
    samples = args["normal_sample"] + "\t" + samples
    libraries = args["normal_library"] + "\t" + libraries

ShellWrapper(snakemake).run(
    r"""
bcftools view \
    --no-update \
    --samples-file <(echo "{libraries}" | tr '\t' '\n') \
    {snakemake.input.annotated} \
    | bcftools reheader \
        --samples <(echo "{samples}" | tr '\t' '\n') \
    | bcftools norm \
        --multiallelics -any \
        --output-type z --output {snakemake.output.vcf} --write-index=tbi
"""
)
