# -*- coding: utf-8 -*-
"""Wrapper for cnvkit.py fix"""

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

# Fix requires empty antitarget file in WGS & Panel modes
if snakemake.input.get("antitarget", None) is None:
    antitarget = "$TMPDIR/antitarget.bed"
    create_dummy_antitarget = f"touch {antitarget} ; "
else:
    antitarget = snakemake.input.antitarget
    create_dummy_antitarget = ""

cmd = r"""
cnvkit.py fix \
    -o {snakemake.output.ratios} \
    {cluster} --sample-id {args[sample-id]} \
    {no_gc} {no_edge} {no_rmask} \
    {snakemake.input.target} {antitarget} {snakemake.input.reference}
""".format(
    snakemake=snakemake,
    args=args,
    antitarget=antitarget,
    cluster="--cluster" if args.get("cluster", False) else "",
    no_gc="--no-gc" if args.get("no-gc", False) else "",
    no_edge="--no-edge" if args.get("no-edge", False) else "",
    no_rmask="--no-rmask" if args.get("no-rmask", False) else "",
)

CnvkitWrapper(snakemake, create_dummy_antitarget + cmd).run()
