# -*- coding: utf-8 -*-
"""Wrapper vor cnvkit.py segment"""

from snappy_wrappers.snappy_wrapper import ShellWrapper

__author__ = "Manuel Holtgrewe"
__email__ = "manuel.holtgrewe@bih-charite.de"

args = getattr(snakemake.params, "args", {})

method = args["method"]
if method == "cbs" and args["smooth_cbs"]:
    method += " --smooth-cbs"

if float(args["threshold"]) > 0:
    threshold = " --threshold " + str(args["threshold"])
else:
    threshold = ""

ShellWrapper(snakemake).run(
    r"""
set -x

# -----------------------------------------------------------------------------

cnvkit.py segment \
    --output {snakemake.output.segments} \
    --method {method} \
    $(if [[ "{args[drop_low_coverage]}" = "True" ]]; then \
        echo --drop-low-coverage
    fi) \
    {threshold} \
    --drop-outliers {args[drop_outliers]} \
    {snakemake.input}
"""
)
