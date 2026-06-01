# -*- coding: utf-8 -*-
"""Wrapper for cnvkit.py antitarget"""

from snappy_wrappers.snappy_wrapper import ShellWrapper

__author__ = "Manuel Holtgrewe"
__email__ = "manuel.holtgrewe@bih-charite.de"

args = getattr(snakemake.params, "args", {})

target = getattr(snakemake.input, "target", "")
if access := getattr(snakemake.input, "access", ""):
    access = f"--access {access}"

# Avoid compariing with None
args["min_size"] = args["min_size"] if args.get("min_size", None) else 0
args["avg_size"] = args["avg_size"] if args.get("avg_size", None) else 0

ShellWrapper(snakemake).run(
    r"""
set -x

# -----------------------------------------------------------------------------

if [[ -n "{target}" ]]
then
    cnvkit.py antitarget \
        --output {snakemake.output.antitarget} \
        {access} \
        $(if [[ {args[avg_size]} -gt 0 ]]; then \
            echo --avg-size {args[avg_size]}
        fi) \
        $(if [[ {args[min_size]} -gt 0 ]]; then \
            echo --min-size {args[min_size]}
        fi) \
        {target}
else
    touch {snakemake.output.antitarget}
fi
"""
)
