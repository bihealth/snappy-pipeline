# -*- coding: utf-8 -*-
"""Wrapper for cnvkit.py coverage"""

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
cnvkit.py coverage --processes {snakemake.resources._cores} \
    -o {snakemake.output.coverage} \
    --fasta {args['reference']} \
    --min-mapq {args[min_mapq]} {count} \
    {args['bam']} {args['intervals']}
""".format(
    snakemake=snakemake,
    args=args,
    count="--count" if args.get("count", False) else "",
)

CnvkitWrapper(snakemake, cmd).run()
