# -*- coding: utf-8 -*-
"""Snakemake wrapper for running ``snappy_wrappers.tools.vcf_filter_denovo.main()``.

isort:skip_file
"""

import collections
from typing import TYPE_CHECKING

from snappy_wrappers.snappy_wrapper import ShellWrapper

if TYPE_CHECKING:
    from snakemake.iocontainers import snakemake

import snappy_wrappers.tools.vcf_filter_denovo

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

args = getattr(snakemake.params, "args", {})

# Build arguments ==================================================================================

besenbacher = args["params_besenbacher"]

# Define arguments in a dictionary, will convert to namedtuple below.
args_all = {
    "verbose": True,
    "index_name": args["index_library"],
    "input_vcf": snakemake.input.vcf,
    "input_ped": snakemake.input.ped,
    "output_vcf": snakemake.output.vcf,
    "offspring_bam": snakemake.input.bam,
    "regions": [],  # empty => process all
    "skip_invalid": False,
    "exclusive_neighborhood": 1000,
    "mnv_neighborhood": 20,
    "use_phase_info": True,
    "haplotype_window": 100000,
    "phase_paternal_first": True,
}
# Bulk-add besenbacher parameters
args_all.update(args["params_besenbacher"])

args_t = collections.namedtuple("Arguments", args_all.keys())(**args_all)

# Execute filtration ===============================================================================

snappy_wrappers.tools.vcf_filter_denovo.run(args_t)

# Postprocess result ===============================================================================

ShellWrapper(snakemake).run(
    r"""
# Build tabix index
tabix -f {snakemake.output.vcf}
"""
)
