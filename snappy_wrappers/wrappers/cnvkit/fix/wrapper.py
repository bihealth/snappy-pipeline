# -*- coding: utf-8 -*-
"""Wrapper for cnvkit.py fix"""

from snappy_wrappers.snappy_wrapper import ShellWrapper

__author__ = "Manuel Holtgrewe"
__email__ = "manuel.holtgrewe@bih-charite.de"

args = getattr(snakemake.params, "args", {})

gender = " --gender {}".format(args["gender"]) if args.get("gender", None) else ""
male = " --male-reference" if args.get("male_reference", False) else ""
no_gc = " --no-gc" if not args["gc_correction"] else ""
no_edge = " --no-edge" if not args["edge_correction"] else ""
no_rmask = " --no-rmask" if not args["rmask_correction"] else ""

ShellWrapper(snakemake).run(
    r"""
set -x

# -----------------------------------------------------------------------------

cnvkit.py fix \
    --output {snakemake.output.ratios} \
    {gender} {male} {no_gc} {no_edge} {no_rmask} \
    {snakemake.input.target} \
    {snakemake.input.antitarget} \
    {snakemake.input.ref}
"""
)
