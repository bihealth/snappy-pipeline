# -*- coding: utf-8 -*-
"""Wrapper for cnvkit.py access"""

from snappy_wrappers.snappy_wrapper import ShellWrapper

__author__ = "Manuel Holtgrewe"
__email__ = "manuel.holtgrewe@bih-charite.de"

args = getattr(snakemake.params, "args", {})
config = args.get("config", {})

exclude = " --exclude " + " -x ".join(config["exclude"]) if config["exclude"] else ""

ShellWrapper(snakemake).run(
    r"""
set -x

# -----------------------------------------------------------------------------

cnvkit.py access \
    -o {snakemake.output.access} \
    $(if [[ {config[min_gap_size]} -gt 0 ]]; then \
        echo --min-gap-size {config[min_gap_size]}
    fi) \
    {exclude} \
    {snakemake.input.reference}
"""
)
