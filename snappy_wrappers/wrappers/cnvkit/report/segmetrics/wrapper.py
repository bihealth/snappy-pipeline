# -*- coding: utf-8 -*-
"""Wrapper for cnvkit.py segmetrics"""

import os
import sys

# The following is required for being able to import snappy_wrappers modules
# inside wrappers.  These run in an "inner" snakemake process which uses its
# own conda environment which cannot see the snappy_pipeline installation.
base_dir = os.path.normpath(os.path.join(os.path.dirname(__file__), "..", "..", "..", ".."))
sys.path.insert(0, base_dir)

from snappy_wrappers.wrappers.cnvkit.cnvkit_wrapper import CnvkitWrapper

args = snakemake.params.get("args", {})

cmd = r"""
cnvkit.py segmetrics \
    -o {snakemake.output.report} \
    --segment {args['segments']} \
    --alpha {args['alpha']} --bootstrap {args['bootstrap']} {smooth_bootstrap} \
    {drop_low_coverage} \
    {stats} \
    {args[ratios]}
""".format(
    snakemake=snakemake,
    args=args,
    drop_low_coverage="--drop-low-coverage" if args.get("drop-low-coverage", False) else "",
    smooth_bootstrap="--smooth-bootstrap" if args.get("smooth-bootstrap", False) else "",
    stats=" ".join([f"--{stat}" for stat in args["stats"]]),
)

CnvkitWrapper(snakemake, cmd).run()
