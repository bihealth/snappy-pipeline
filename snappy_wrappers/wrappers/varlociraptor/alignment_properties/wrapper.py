# -*- coding: utf-8 -*-
"""Estimate the alignment properties of a library for ``varlociraptor preprocess``."""

from typing import TYPE_CHECKING

from snappy_wrappers.snappy_wrapper import ShellWrapper

if TYPE_CHECKING:
    from snakemake.iocontainers import snakemake

ShellWrapper(snakemake).run(
    r"""
varlociraptor estimate alignment-properties \
    {snakemake.input.reference} \
    --bams {snakemake.input.bam} \
    > {snakemake.output.json}
"""
)
