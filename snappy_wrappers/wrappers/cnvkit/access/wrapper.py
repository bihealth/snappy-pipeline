# -*- coding: utf-8 -*-
"""Wrapper for cnvkit.py access"""

import os
import sys

# The following is required for being able to import snappy_wrappers modules
# inside wrappers.  These run in an "inner" snakemake process which uses its
# own conda environment which cannot see the snappy_pipeline installation.
base_dir = os.path.normpath(os.path.join(os.path.dirname(__file__), "..", "..", "..", ".."))
sys.path.insert(0, base_dir)

from snappy_wrappers.tools.genome_windows import ignore_chroms
from snappy_wrappers.wrappers.cnvkit.cnvkit_wrapper import CnvkitWrapper

__author__ = "Eric Blanc"
__email__ = "eric.blanc@bih-charite.de"

args = snakemake.params.get("args", {})

prefix = ""

# Add the "ignore_chrom" contents to the excluded regions
if len(args.get("ignore_chroms", [])) > 0:
    ignored_contigs = ignore_chroms(args["reference"], args["ignore_chroms"], return_ignored=True)
    lines = ["cat << __EOF > $TMPDIR/ignore_chroms.bed"]
    for (contig_name, contig_length) in ignored_contigs:
        lines.append(f"{contig_name}\t0\t{contig_length}")
    lines.append("__EOF")
    prefix = "\n".join(lines) + "\n"
    args["exclude"].append("$TMPDIR/ignore_chroms.bed")

cmd = r"""
cnvkit.py access \
    -o {snakemake.output.access} \
    {min_gap_size} {exclude} \
    {args[reference]}
""".format(
    snakemake=snakemake,
    args=args,
    min_gap_size=f"--min-gap-size {args['min-gap-size']}" if args.get("min-gap-size", None) is not None else "",
    exclude=" ".join([f"--exclude {x}" for x in args["exclude"]]),
)

CnvkitWrapper(snakemake, prefix + cmd).run()
