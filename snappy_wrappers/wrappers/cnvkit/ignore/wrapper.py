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

ignored_contigs = ignore_chroms(snakemake.input.reference, args["ignore_chroms"], return_ignored=True)
lines = []
for (contig_name, contig_length) in ignored_contigs:
    lines.append(f"{contig_name}\t0\t{contig_length}")
lines = "\n".join(lines)

cmd = f"""
cat << __EOF > {snakemake.output.ignore_chroms}
{lines}
__EOF
"""

CnvkitWrapper(snakemake, cmd).run()
