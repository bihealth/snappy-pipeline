# -*- coding: utf-8 -*-
"""Snakemake wrapper for running ``snappy_wrappers.tools.vcf_filter_denovo.main()``.

isort:skip_file
"""

import snappy_wrappers.tools.vcf_filter_denovo
from snakemake.shell import shell
import collections
import os
import sys


__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

# Build arguments ==================================================================================

besenbacher = snakemake.config["step_config"]["variant_denovo_filtration"]["params_besenbacher"]

# Define arguments in a dictionary, will convert to namedtuple below.
args = {
    "verbose": True,
    "index_name": snakemake.wildcards.index_library,
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
args.update(snakemake.config["step_config"]["variant_denovo_filtration"]["params_besenbacher"])

args_t = collections.namedtuple("Arguments", args.keys())(**args)

# Execute filtration ===============================================================================

snappy_wrappers.tools.vcf_filter_denovo.run(args_t)

# Postprocess result ===============================================================================

shell(
    r"""
# Build tabix index
tabix -f {snakemake.output.vcf}

# Compute MD5 sums
pushd $(dirname {snakemake.output.vcf})
md5sum $(basename {snakemake.output.vcf}) > $(basename {snakemake.output.vcf}).md5
md5sum $(basename {snakemake.output.tbi}) > $(basename {snakemake.output.tbi}).md5
"""
)
