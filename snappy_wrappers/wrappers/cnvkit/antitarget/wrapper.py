# -*- coding: utf-8 -*-
"""Wrapper for cnvkit.py antitarget"""

from snappy_wrappers.wrappers.cnvkit.cnvkit_wrapper import CnvkitWrapper

__author__ = "Eric Blanc"
__email__ = "eric.blanc@bih-charite.de"

cmd = r"""
cnvkit.py antitarget \
    -o {snakemake.output.region} \
    --avg-size {snakemake.params.avg_size} --min-size {snakemake.params.min_size} \
    {access} \
    {snakemake.input.target}
""".format(
    snakemake=snakemake,
    access=f"--access {snakemake.params.access}" if snakemake.params.access else "",
)

CnvkitWrapper(snakemake, cmd).run()
