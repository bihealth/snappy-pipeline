# -*- coding: utf-8 -*-
"""
Wrapper for annotating germline variants using dbSNP, gnomAD, ... or with another vcf

It is meant to be used in conjunction with other bcftools commands, such as call & filter

Mandatory snakemake.input: vcf
Optional snakemake.input: annotations, columns_file, header_lines, regions_file, samples_file


Mandatory snakemake.params.args: extra_args
Optional snakemake.params.args: index

Mandatory snakemake.output: vcf
Optional snakemake.output:
"""

__author__ = "Eric Blanc <eric.blanc@bih-charite.de>"

import os
import sys

# The following is required for being able to import snappy_wrappers modules
# inside wrappers. When the wrappers have their own python environment, messing
# with the path is necessary.
base_dir = os.path.normpath(os.path.join(os.path.dirname(__file__), "..", "..", "..", "..", ".."))
sys.path.insert(0, base_dir)

from snappy_wrappers.snappy_wrapper import BcftoolsWrapper

wrapper = BcftoolsWrapper(snakemake)
cmd = wrapper.get_command(tool="annotate")
wrapper.run(cmd)
