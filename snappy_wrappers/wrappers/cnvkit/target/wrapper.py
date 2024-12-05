# -*- coding: utf-8 -*-
"""Wrapper for cnvkit.py target"""

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
cnvkit.py target \
    -o {snakemake.output.target} \
    {avg_size} {split} {annotate} {short_names} \
    {snakemake.input.interval}
""".format(
    snakemake=snakemake,
    args=args,
    avg_size=f"--avg-size {args['avg-size']}" if args.get("avg-size", None) is not None else "",
    split=f"--split" if args.get("split", False) else "",
    annotate=f"--annotate {snakemake.input.annotate}" if snakemake.input.get("annotate", None) is not None else "",
    short_names="--short-names" if args.get("short-names", False) else "",
)

CnvkitWrapper(snakemake, cmd).run()
