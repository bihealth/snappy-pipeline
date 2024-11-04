# -*- coding: utf-8 -*-
"""Wrapper for cnvkit.py fix"""

from snappy_wrappers.wrappers.cnvkit.cnvkit_wrapper import CnvkitWrapper

__author__ = "Eric Blanc"
__email__ = "eric.blanc@bih-charite.de"

cmd = r"""
cnvkit.py fix \
    -o {snakemake.output.coverage} \
    {cluster} {snakemake.params.sample_id} \
    {no_gc} {no_edge} {no_rmask} \
    {snakemake.input.target} {antitarget} {snakemake.input.reference}
""".format(
    snakemake=snakemake,
    cluster="--cluster" if snakemake.params.cluster else "",
    no_gc="--no-gc" if snakemake.params.no_gc else "",
    no_edge="--no-edge" if snakemake.params.no_edge else "",
    no_rmask="--no-rmask" if snakemake.params.no_rmask else "",
    antitarget=f"{snakemake.input.antitarget}" if "antitarget" in snakemake.input else "",
)

CnvkitWrapper(snakemake, cmd).run()
