# -*- coding: utf-8 -*-
"""Wrapper vor cnvkit.py segment"""

import os
import sys

# The following is required for being able to import snappy_wrappers modules
# inside wrappers.  These run in an "inner" snakemake process which uses its
# own conda environment which cannot see the snappy_pipeline installation.
base_dir = os.path.normpath(os.path.join(os.path.dirname(__file__), "..", "..", "..", ".."))
sys.path.insert(0, base_dir)

from snappy_wrappers.wrappers.cnvkit.cnvkit_wrapper import CnvkitWrapper

args = snakemake.params.get("args", {})

if snakemake.input.get("variants", None) is not None:
    variants = r"""
        ---vcf {snakemake.input.variants} \
        --sample-id {args[sample-id]} --normal-id {args[normal-id]} \
        --min-variant-depth {args[min-variant-depth]} {zygocity_freq}
    """.format(
        snakemake=snakemake,
        args=args,
        zygocity_freq=f"--zygocity_freq {args['zygocity-freq']}" if "zygocity-freq" in args else ""
    )
else:
    variants = ""

cmd = r"""
cnvkit.py segment --processes {snakemake.resources._cores} \
    -o {snakemake.output.segments} --dataframe {snakemake.output.dataframe} \
    --method {args[method]} --threshold {args[threshold]} {smooth_cbs} \
    {drop_low_coverage} --drop-outliers {args[drop-outliers]} \
    {variants} \
    {snakemake.input.ratios}
""".format(
    snakemake=snakemake,
    args=args,
    variants=variants,
    smooth_cbs="--smooth-cbs" if args.get("smooth-cbs", False) else "",
    drop_low_coverage="--drop-low-coverage" if args.get("drop-low-coverage", False) else "",
)

CnvkitWrapper(snakemake, cmd).run()
