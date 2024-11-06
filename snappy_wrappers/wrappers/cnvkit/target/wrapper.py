# -*- coding: utf-8 -*-
"""Wrapper for cnvkit.py target"""

import os
import re
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

# WGS: targets are all accessible regions, WES: targets are baits
interval = snakemake.input.access if snakemake.input.get("access", None) else args["target"]

if snakemake.input.get("avg_size", "") != "":
    pattern = re.compile("^Target:[ \t]+([+-]?(\d+(\.\d*)?|\.\d+)([EeDd][+-]?[0-9]+)?)[ \t]+([+-]?(\d+(\.\d*)?|\.\d+)([EeDd][+-]?[0-9]+)?)$")
    with open(snakemake.input.avg_size) as f:
        for line in f:
            m = pattern.match(line)
            if m:
                avg_size = int(float(m.groups()[4]))
                break

else:
    avg_size = args["avg_size"]

cmd = r"""
cnvkit.py target \
    -o {snakemake.output.target} \
    {avg_size} {split} {annotate} \
    {interval}
""".format(
    snakemake=snakemake,
    args=args,
    interval=interval,
    avg_size=f"--avg-size {avg_size}",
    split=f"--split" if "split" in args and args["split"] else "",
    annotate=f"--annotate {args['annotate']}" if "annotate" in args else "",
)

CnvkitWrapper(snakemake, cmd).run()
