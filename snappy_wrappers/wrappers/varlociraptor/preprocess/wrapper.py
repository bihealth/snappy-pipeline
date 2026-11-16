# -*- coding: utf-8 -*-
"""Observations of a library at the candidates of one genome window.

A candidate belongs to the window of its position, so the windows of an entity partition its
candidates and every library of the entity sees the same candidates in a window.
"""

from typing import TYPE_CHECKING

from snappy_wrappers.snappy_wrapper import ShellWrapper

if TYPE_CHECKING:
    from snakemake.iocontainers import snakemake

args = snakemake.params.args

ShellWrapper(snakemake).run(
    r"""
bcftools view \
    --regions-file {snakemake.input.region} \
    --regions-overlap pos \
    --output-type b \
    --output $TMPDIR/candidates.bcf \
    {snakemake.input.candidates}

varlociraptor preprocess variants \
    {snakemake.input.reference} \
    --candidates $TMPDIR/candidates.bcf \
    --bam {snakemake.input.bam} \
    --alignment-properties {snakemake.input.alignment_properties} \
    --max-depth {args[max_depth]} \
    {extra_arguments} \
    --output {snakemake.output.bcf}
"""
)
