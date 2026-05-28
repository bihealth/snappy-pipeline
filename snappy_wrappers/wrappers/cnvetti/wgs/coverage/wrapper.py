# -*- coding: utf-8 -*-
"""Wrapper for running CNVetti WGS coverage step."""

from snappy_wrappers.snappy_wrapper import ShellWrapper

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

args = getattr(snakemake.params, "args", {})

ShellWrapper(snakemake).run(
    r"""
set -x

REF={args[reference]}

cnvetti cmd coverage \
    -vvv \
    --considered-regions GenomeWide \
    --count-kind {args[count_kind]} \
    --window-length {args[window_length]} \
    --reference $REF \
    --output $TMPDIR/cov.bcf \
    --input {snakemake.input.bam}

cnvetti cmd normalize \
    -vvv \
    --normalization {args[normalization]} \
    --input $TMPDIR/cov.bcf \
    --output {snakemake.output.bcf}
"""
)
