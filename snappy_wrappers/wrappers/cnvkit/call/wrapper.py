# -*- coding: utf-8 -*-
"""Wrapper vor cnvkit.py call"""

import os
import re
import sys

# The following is required for being able to import snappy_wrappers modules
# inside wrappers.  These run in an "inner" snakemake process which uses its
# own conda environment which cannot see the snappy_pipeline installation.
base_dir = os.path.normpath(os.path.join(os.path.dirname(__file__), "..", "..", "..", ".."))
sys.path.insert(0, base_dir)

from snappy_wrappers.wrappers.cnvkit.cnvkit_wrapper import CnvkitWrapper

args = snakemake.params.get("args", {})

PATTERN = re.compile("^(Purity|Ploidy): +([+-]?([0-9]+(\.[0-9]*)?|\.[0-9]+)([EeDd][+-]?[0-9]+)?) *$")


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
cnvkit.py call \
    -o {snakemake.output.calls} \
    --method {args[method]} {thresholds} \
    {filter} \
    {center} {center_at} {drop_low_coverage} {sample_sex} {male_reference} \
    {purity} {ploidy} \
    {variants} \
    {args[segments]}
""".format(
    snakemake=snakemake,
    args=args,
    variants=variants,
    purity=f"--purity {args['purity']}" if args.get("purity", None) is not None else "",
    ploidy=f"--ploidy {args['ploidy']}" if args.get("ploidy", None) is not None else "",
    thresholds="--thresholds={}".format(",".join(map(str, args["thresholds"]))) if len(args.get("thresholds", [])) > 0 else "",
    filter=f"--filter {args['filter']}" if args.get("filter", None) is not None else "",
    center=f"--center {args['center']}" if args.get("center", None) is not None else "",
    center_at=f"--center-at {args['center-at']}" if args.get("center-at", None)  is not None else "",
    drop_low_coverage="--drop-low-coverage" if args.get("drop-low-coverage", False) else "",
    sample_sex=f"--sample-sex {args['sample-sex']}" if args.get("sample-sex", None) is not None else "",
    male_reference=f"--male-reference" if args.get("male-reference", False) else "",
)

CnvkitWrapper(snakemake, cmd).run()
