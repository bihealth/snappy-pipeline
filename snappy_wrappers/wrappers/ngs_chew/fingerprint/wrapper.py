# -*- coding: utf-8 -*-
"""Wrapper for running ``ngs-chew fingerprint``."""

from snakemake.shell import shell

from snappy_wrappers.snappy_wrapper import ShellWrapper

__author__ = "Manuel Holtgrewe"
__email__ = "manuel.holtgrewe@bih-charite.de"

args = getattr(snakemake.params, "args", {})

ShellWrapper(snakemake).run(
    r"""
set -x

mkdir -p $TMPDIR/{{out,sorted,sort.tmp}}

ngs-chew fingerprint \
    --reference {args[reference]} \
    --output-aafs \
    --output-fingerprint {snakemake.output.npz} \
    --input-bam {snakemake.input.bam}
"""
)
