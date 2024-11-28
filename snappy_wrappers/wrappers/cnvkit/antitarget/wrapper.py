# -*- coding: utf-8 -*-
"""Wrapper for cnvkit.py antitarget"""

import os
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
cnvkit.py antitarget \
    -o {snakemake.output.antitarget} \
    --avg-size {args[avg-size]} {min_size} \
    --access {args[access]} \
    {args[target]}
""".format(
    snakemake=snakemake,
    args=args,
    min_size=f"--min-size {args['min-size']}" if args.get("min-size") is not None else "",
)

CnvkitWrapper(snakemake, cmd).run()
