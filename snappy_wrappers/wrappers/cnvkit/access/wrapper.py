# -*- coding: utf-8 -*-
"""Wrapper for cnvkit.py access"""

from snappy_wrappers.snappy_wrapper import ShellWrapper

__author__ = "Manuel Holtgrewe"
__email__ = "manuel.holtgrewe@bih-charite.de"

args = getattr(snakemake.params, "args", {})
min_gap_size = args["min_gap_size"]

exclude_files = snakemake.input.get("exclude", [])
exclude = " --exclude " + " -x ".join(exclude_files) if exclude_files else ""

ShellWrapper(snakemake).run(
    r"""
set -x

# -----------------------------------------------------------------------------

cnvkit.py access \
    -o {snakemake.output.access} \
    $(if [[ {min_gap_size} -gt 0 ]]; then \
        echo --min-gap-size {min_gap_size}
    fi) \
    {exclude} \
    {snakemake.input.reference}
"""
)
