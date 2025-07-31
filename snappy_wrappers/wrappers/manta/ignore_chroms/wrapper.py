# -*- coding: utf-8 -*-
"""Wrapper for running Manta in somatic variant calling mode on WGS data"""

import os

from snappy_wrappers.snappy_wrapper import ShellWrapper
from snappy_wrappers.tools.genome_windows import ignore_chroms

__author__ = "Eric Blanc"
__email__ = "eric.blanc@bih-charite.de"

args = snakemake.params.get("args", {})

bed = str(snakemake.output.callRegions).replace(".gz", "")

with open(bed, "wt") as f:
    for contig, length in ignore_chroms(str(snakemake.input.reference), set(args["ignore_chroms"]), return_ignored=False):
        print(f"{contig}\t0\t{length}", file=f)

cmd=f"bgzip {bed} ; tabix {snakemake.output.callRegions}"

ShellWrapper(snakemake).run(cmd)
