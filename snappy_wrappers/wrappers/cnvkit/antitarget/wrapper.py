# -*- coding: utf-8 -*-
"""Wrapper for cnvkit.py antitarget"""

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

if snakemake.input.get("target", "") != "":
    cmd = r"""
    cnvkit.py antitarget \
        -o {snakemake.output.antitarget} \
        --avg-size {args['avg_size']} --min-size {args['min_size']} \
        {access} \
        {snakemake.input.target}
    """.format(
        snakemake=snakemake,
        args=args,
        access=f"--access {args['access']}" if "access" in args else "",
    )
else:
    cmd = f"touch {snakemake.output.antitarget}"

CnvkitWrapper(snakemake, cmd).run()
