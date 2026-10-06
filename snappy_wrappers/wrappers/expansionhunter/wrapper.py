# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for ExpansionHunter: Snakemake wrapper.py"""

import os

from snakemake import shell

shell.executable("/bin/bash")

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


shell(
    r"""
set -x

# TODO: remove this again, is for fail early
# Additional logging for transparency & reproducibility
# Logging: Save a copy this wrapper (with the pickle details in the header)
cp {this_file} $(dirname {snakemake.log})/wrapper_expansionhunter.py

# Also pipe stderr to log file
if [[ -n "{snakemake.log}" ]]; then
    if [[ "$(set +e; tty; set -e)" != "" ]]; then
        rm -f "{snakemake.log}" && mkdir -p $(dirname {snakemake.log})
        exec 2> >(tee -a "{snakemake.log}" >&2)
    else
        rm -f "{snakemake.log}" && mkdir -p $(dirname {snakemake.log})
        echo "No tty, logging disabled" >"{snakemake.log}"
    fi
fi

# Create out dir
mkdir -p $(dirname {snakemake.output.json})

# Call tool
ExpansionHunter --reads {snakemake.input.bam} \
        --reference {snakemake.input.reference} \
        --variant-catalog {snakemake.input.repeat_catalog} \
        --output-prefix {prefix} {sex_argument}
"""
)

# Compute MD5 sums of log and vcf.
shell(
    r"""
md5sum {snakemake.log} > {snakemake.log}.md5
md5sum {snakemake.output.vcf} > {snakemake.output.vcf_md5}
"""
)
