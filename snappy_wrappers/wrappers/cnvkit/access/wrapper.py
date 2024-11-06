# -*- coding: utf-8 -*-
"""Wrapper for cnvkit.py access"""

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
cnvkit.py access \
    -o {snakemake.output.access} \
    --min-gap-size {args[min_gap_size]} \
    {exclude} \
    {args[reference]}
""".format(
    snakemake=snakemake,
    args=args,
    exclude=" ".join([f"--exclude {x}" for x in args["exclude"]]) if "exclude" in args else "",
)

CnvkitWrapper(snakemake, cmd).run()
