# -*- coding: utf-8 -*-
"""Wrapper for cnvkit.py scatter"""

import os
import re
import sys

# The following is required for being able to import snappy_wrappers modules
# inside wrappers.  These run in an "inner" snakemake process which uses its
# own conda environment which cannot see the snappy_pipeline installation.
base_dir = os.path.normpath(os.path.join(os.path.dirname(__file__), "..", "..", "..", "..", ".."))
sys.path.insert(0, base_dir)

from snappy_wrappers.wrappers.cnvkit.cnvkit_wrapper import CnvkitWrapper

args = snakemake.params.get("args", {})

if "variants" in args:
    variants = r"""
        ---vcf {args[variants]} \
        --sample-id {args[sample-id]} --normal-id {args[normal-id]} \
        --min-variant-depth {args[min-variant-depth]} {zygocity_freq}
    """.format(
        snakemake=snakemake,
        args=args,
        zygocity_freq=f"--zygocity-freq {args['zygocity-freq']}" if args.get("zygocity-freq", None) is not None else "",
    )
else:
    variants = ""

cmd = r"""
cnvkit.py scatter \
    -o {snakemake.output.plot} \
    --segment {args[segments]} \
    {chromosome} {gene} {range_list} \
    --width {args[width]} \
    --antitarget-marker {args[antitarget-marker]} --segment-color {args[segment-color]} \
    {by_bin} {trend} --title "{args[title]}" \
    {y_min} {y_max} {fig_size} \
    {variants} \
    {args[ratios]}
""".format(
    snakemake=snakemake,
    args=args,
    variants=variants,
    chromosome=f"--chromosome {args['chromosome']}" if args.get("chromosome", None) is not None else "",
    gene=f"--gene {args['gene']}" if args.get("gene", None) is not None else "",
    range_list=f"--range-list {args['range-list']}" if args.get("range-list", None) is not None else "",
    by_bin="--by-bin" if args.get("by-bin", False) else "",
    trend="--trend" if args.get("trend", False) else "",
    y_min=f"--y-min {args['y-min']}" if args.get("y-min", None) is not None else "",
    y_max=f"--y-max {args['y-max']}" if args.get("y-max", None) is not None else "",
    fig_size="--fig-size {}".format(" ".join(map(str, args['fig-size']))) if args.get("fig-size", None) is not None else "",
)

CnvkitWrapper(snakemake, cmd).run()
