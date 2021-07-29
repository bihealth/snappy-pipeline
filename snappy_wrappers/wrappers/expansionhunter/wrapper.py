# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for ExpansionHunter: Snakemake wrapper.py
"""
import os

from snakemake import shell

shell.executable("/bin/bash")

this_file = __file__

# Define prefix based on json output
prefix = snakemake.output.json
prefix = prefix.replace(".json", "")
prefix = os.path.join(os.getcwd(), prefix)


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

ExpansionHunter --reads {snakemake.input} \
        --reference {snakemake.config[static_data_config][reference][path]} \
        --variant-catalog {snakemake.config[step_config][repeat_expansion][repeat_catalog]} \
        --output-prefix {prefix}


"""
)

# Compute MD5 sums of logs.
shell(
    r"""
md5sum {snakemake.log} >{snakemake.log}.md5
"""
)
