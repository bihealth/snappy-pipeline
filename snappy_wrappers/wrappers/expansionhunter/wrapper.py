# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for ExpansionHunter: Snakemake wrapper.py"""

import os
from typing import TYPE_CHECKING

from snappy_wrappers.snappy_wrapper import ShellWrapper

if TYPE_CHECKING:
    from snakemake.iocontainers import snakemake

this_file = __file__

args = getattr(snakemake.params, "args", {})

# Define prefix based on json output
prefix = snakemake.output.json
prefix = prefix.replace(".json", "")
prefix = os.path.join(os.getcwd(), prefix)

# Define argument sex if any (otherwise: female [default])
sex_argument = ""
valid_sex_list = ["female", "male"]
if args["sex"] in valid_sex_list:
    sex_argument = "--sex " + args["sex"]


ShellWrapper(snakemake).run(
    r"""
# Additional logging for transparency & reproducibility
# Logging: Save a copy this wrapper (with the pickle details in the header)
cp {this_file} $(dirname {snakemake.log})/wrapper_expansionhunter.py

# Create out dir
mkdir -p $(dirname {snakemake.output.json})

# Call tool
ExpansionHunter --reads {snakemake.input.bam} \
        --reference {snakemake.input.reference} \
        --variant-catalog {snakemake.input.repeat_catalog} \
        --output-prefix {prefix} {sex_argument}
"""
)
