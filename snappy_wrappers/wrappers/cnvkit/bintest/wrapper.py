# -*- coding: utf-8 -*-
"""Wrapper for cnvkit.py bintest"""

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
cnvkit.py bintest \
    -o {snakemake.output.tests} \
    --segment {args[segments]} \
    --alpha {args[alpha]} {target} \
    {args[ratios]}
""".format(
    snakemake=snakemake,
    args=args,
    target=f"--target" if args.get("target", False) else "",
)

CnvkitWrapper(snakemake, cmd).run()
