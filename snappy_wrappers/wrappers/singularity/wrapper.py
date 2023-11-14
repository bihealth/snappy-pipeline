# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for preparing exome kit intervals for PureCN
"""

from snakemake import shell

__author__ = "Eric Blanc <eric.blanc@bih-charite.de>"

container = None
if "container" in snakemake.params.keys() and snakemake.params["container"]:
    container = snakemake.params["container"]
else:
    step = snakemake.config["pipeline_step"]["name"]
    config = snakemake.config["step_config"][step]
    if "container" in config.keys() and config["container"]:
        container = config["container"]
assert container, "Missing or illegal container image address"

shell.executable("/bin/bash")

shell(
    r"""
set -x

# Also pipe everything to log file
if [[ -n "{snakemake.log.log}" ]]; then
    if [[ "$(set +e; tty; set -e)" != "" ]]; then
        rm -f "{snakemake.log.log}" && mkdir -p $(dirname {snakemake.log.log})
        exec &> >(tee -a "{snakemake.log.log}" >&2)
    else
        rm -f "{snakemake.log.log}" && mkdir -p $(dirname {snakemake.log.log})
        echo "No tty, logging disabled" >"{snakemake.log.log}"
    fi
fi

# Write out information about conda installation.
conda list >{snakemake.log.conda_list}
conda info >{snakemake.log.conda_info}
md5sum {snakemake.log.conda_list} >{snakemake.log.conda_list_md5}
md5sum {snakemake.log.conda_info} >{snakemake.log.conda_info_md5}

apptainer pull --name {snakemake.output.container} {container}
"""
)

# Compute MD5 sums of logs.
shell(
    r"""
md5sum {snakemake.log.log} >{snakemake.log.log_md5}
"""
)
