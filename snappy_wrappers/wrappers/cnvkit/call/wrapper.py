# -*- coding: utf-8 -*-
"""Wrapper vor cnvkit.py call"""

from snappy_wrappers.snappy_wrapper import ShellWrapper

__author__ = "Manuel Holtgrewe"
__email__ = "manuel.holtgrewe@bih-charite.de"

args = getattr(snakemake.params, "args", {})

if center := args.get("center", None):
    if center in ("mean", "median", "mode", "biweight"):
        center = " --center " + center
    else:
        center = " --center-at " + str(center)
else:
    center = ""

gender = " --gender {}".format(args["gender"]) if args.get("gender", None) else ""
male = " --male-reference" if args.get("male_reference", False) else ""
purity = " --purity {}".format(args["purity"]) if args.get("purity", 0.0) > 0.0 else ""

ShellWrapper(snakemake).run(
    r"""
set -x

# -----------------------------------------------------------------------------

cnvkit.py call \
    --output {snakemake.output.calls} \
    --method {args[method]} \
    --thresholds={args[thresholds]} \
    $(if [[ -n "{args[filter]}" ]]; then \
        echo --filter {args[filter]}
    fi) \
    {center} {gender} {male} \
    --ploidy {args[ploidy]} {purity} \
    {snakemake.input}
"""
)
