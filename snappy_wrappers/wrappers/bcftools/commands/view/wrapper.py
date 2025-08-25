# -*- coding: utf-8 -*-
"""
Wrapper for indexing germline variants

It is meant to be used in conjunction with other bcftools commands, such as mpileup & call

Mandatory snakemake.input: vcf
Optional snakemake.input: regions_file, samples_file, targets_file

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
base_dir = os.path.normpath(os.path.dirname(__file__))
while os.path.basename(base_dir) != "snappy_wrappers":
    base_dir = os.path.dirname(base_dir)
sys.path.insert(0, os.path.dirname(base_dir))

from snappy_wrappers.snappy_wrapper import BcftoolsWrapper

wrapper = BcftoolsWrapper(snakemake)
cmd = wrapper.get_command(tool="view")
wrapper.run(cmd)
