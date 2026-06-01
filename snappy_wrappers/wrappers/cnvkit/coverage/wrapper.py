# -*- coding: utf-8 -*-
"""Wrapper for cnvkit.py coverage"""

from snappy_wrappers.snappy_wrapper import ShellWrapper

__author__ = "Manuel Holtgrewe"
__email__ = "manuel.holtgrewe@bih-charite.de"

args = getattr(snakemake.params, "args", {})

ShellWrapper(snakemake).run(
    r"""
set -x

# Function definitions ---------------------------------------------------------

coverage()
{{
    cnvkit.py coverage \
        --fasta {snakemake.input.reference} \
        --min-mapq {args[min_mapq]} \
        --processes {snakemake.threads} \
        {snakemake.input.bam} \
        --output $2 $1
}}

# -----------------------------------------------------------------------------

coverage {snakemake.input.target} {snakemake.output.target}
coverage {snakemake.input.antitarget} {snakemake.output.antitarget}
"""
)
