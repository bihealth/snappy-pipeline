# -*- coding: utf-8 -*-
"""Wrapper for cnvkit.py genemetrics"""

import os
import sys

# The following is required for being able to import snappy_wrappers modules
# inside wrappers.  These run in an "inner" snakemake process which uses its
# own conda environment which cannot see the snappy_pipeline installation.
base_dir = os.path.normpath(os.path.join(os.path.dirname(__file__), "..", "..", "..", "..", ".."))
sys.path.insert(0, base_dir)

from snappy_wrappers.wrappers.cnvkit.cnvkit_wrapper import CnvkitWrapper

args = snakemake.params.get("args", {})

# As segment files can appear multiple times, the `--segment` argument must be the last one
cmd = r"""
cnvkit.py metrics \
    -o {snakemake.output.report} \
    {drop_low_coverage} \
    {snakemake.input.ratios} \
    {segments}
""".format(
    snakemake=snakemake,
    args=args,
    drop_low_coverage="--drop-low-coverage" if args.get("drop-low-coverage", False) else "",
    segments=f"--segments {snakemake.input.segments}" if snakemake.input.get("segments", None) is not None else "",
)

CnvkitWrapper(snakemake, cmd).run()
