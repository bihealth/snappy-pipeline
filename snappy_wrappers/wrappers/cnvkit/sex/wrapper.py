# -*- coding: utf-8 -*-
"""Wrapper for cnvkit.py sex"""

import os
import re
import sys

# The following is required for being able to import snappy_wrappers modules
# inside wrappers.  These run in an "inner" snakemake process which uses its
# own conda environment which cannot see the snappy_pipeline installation.
base_dir = os.path.normpath(os.path.join(os.path.dirname(__file__), "..", "..", "..", ".."))
sys.path.insert(0, base_dir)

from snappy_wrappers.wrappers.cnvkit.cnvkit_wrapper import CnvkitWrapper

__author__ = "Eric Blanc"
__email__ = "eric.blanc@bih-charite.de"

args = snakemake.params.get("args", {})

cmd = r"""
cnvkit.py sex \
    -o {snakemake.output.sex} \
    {diploid_parx_genome} \
    {coverages}
""".format(
    snakemake=snakemake,
    args=args,
    diploid_parx_genome=f"--diploid-parx-genome {args['diploid-parx-genome']}" if args.get('diploid-parx-genome', None) is not None else "",
    coverages=" ".join(args["coverages"]),
)

CnvkitWrapper(snakemake, cmd).run()
