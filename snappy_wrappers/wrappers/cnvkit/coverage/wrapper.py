# -*- coding: utf-8 -*-
"""Wrapper for cnvkit.py coverage"""

from snappy_wrappers.wrappers.cnvkit.cnvkit_wrapper import CnvkitWrapper

__author__ = "Eric Blanc"
__email__ = "eric.blanc@bih-charite.de"

cmd = r"""
cnvkit.py coverage --processes {snakemake.params.processes} \
    -o {snakemake.output.coverage} \
    --fasta {snakemake.params.reference}
    --min-mapq {snakemake.params.min_mapq} {count} \
    {snakemake.input.bam} {snakemake.input.intervals}
""".format(
    snakemake=snakemake,
    count="--count" if snakemake.params.count else "",
)

CnvkitWrapper(snakemake, cmd).run()
