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
exclude = args.get("exclude", [])

# Add the "ignore_chrom" contents to the excluded regions
if snakemake.input.get("ignore_chroms", None) is not None:
    exclude.append(snakemake.input.get("ignore_chroms"))

cmd = r"""
cnvkit.py access \
    -o {snakemake.output.access} \
    {min_gap_size} {exclude} \
    {snakemake.input.reference}
""".format(
    snakemake=snakemake,
    args=args,
    min_gap_size=f"--min-gap-size {args['min-gap-size']}" if args.get("min-gap-size", None) is not None else "",
    exclude=" ".join([f"--exclude {x}" for x in exclude]),
)

CnvkitWrapper(snakemake, cmd).run()
