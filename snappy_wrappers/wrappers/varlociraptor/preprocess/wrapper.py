# -*- coding: utf-8 -*-
"""Observations of a library at one chunk of the candidates of its entity."""

from typing import TYPE_CHECKING

from snappy_wrappers.snappy_wrapper import ShellWrapper

if TYPE_CHECKING:
    from snakemake.iocontainers import snakemake

args = snakemake.params.args

ShellWrapper(snakemake).run(
    r"""
varlociraptor preprocess variants \
    {snakemake.input.reference} \
    --candidates {snakemake.input.candidates} \
    --bam {snakemake.input.bam} \
    --alignment-properties {snakemake.input.alignment_properties} \
    --max-depth {args[max_depth]} \
    {extra_arguments} \
    --output {snakemake.output.bcf}
"""
)
