# -*- coding: utf-8 -*-
"""Wrapper for cnvkit.py autobin (replicating cnvkit batch)"""

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

cmd = r"""
cnvkit.py autobin --method {args[method]} \
    {out_target} {out_antitarget} \
    {access} {target} \
    --bp-per-bin {args[bp-per-bin]} \
    {snakemake.input.bams} \
    > {snakemake.output.result}
""".format(
    snakemake=snakemake,
    args=args,
    out_target=f"--target-output-bed {snakemake.output.target}" if snakemake.output.get("target", "") != "" else "",
    out_antitarget=f"--antitarget-output-bed {snakemake.output.antitarget}" if snakemake.output.get("antitarget", "") != "" else "",
    access=f"--access {snakemake.input.access}" if snakemake.input.get("access", None) is not None else "",
    target=f"--targets {snakemake.input.target}" if snakemake.input.get("target", None) is not None else "",
)

CnvkitWrapper(snakemake, cmd).run()
