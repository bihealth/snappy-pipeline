# -*- coding: utf-8 -*-
"""Wrapper for running CNVetti WGS segment step."""

from snappy_wrappers.snappy_wrapper import ShellWrapper

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

args = getattr(snakemake.params, "args", {})

ShellWrapper(snakemake).run(
    r"""
set -x

cnvetti cmd segment \
    -vvv \
    --segmentation {args[segmentation]} \
    --input {snakemake.input.bcf} \
    --output {snakemake.output.windows_bcf} \
    --output-segments {snakemake.output.segments_bcf}
"""
)
