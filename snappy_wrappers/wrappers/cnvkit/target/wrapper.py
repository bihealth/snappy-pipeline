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

# WGS: targets are all accessible regions, WES: targets are baits
interval = snakemake.input.access if "access" in snakemake.input else snakemake.params.target

if "avg_size" in snakemake.input:
    pattern = re.compile("^Target:[ \t]+([+-]?(\d+(\.\d*)?|\.\d+)([EeDd][+-]?[0-9]+)?)[ \t]+([+-]?(\d+(\.\d*)?|\.\d+)([EeDd][+-]?[0-9]+)?)$")
    with open(snakemake.input.avg_size) as f:
        for line in f:
            m = pattern.match(line)
            if m:
                avg_size = float(m.groups()[4])
                break

else:
    avg_size = snakemake.params.avg_size

cmd = r"""
cnvkit.py target \
    -o {snakemake.output.target} \
    {avg_size} {split} \
    {interval}
""".format(
    snakemake=snakemake,
    interval=interval,
    avg_size=f"--avg-size {avg_size}",
    split=f"--split" if snakemake.params.split else "",
)

CnvkitWrapper(snakemake, cmd).run()
