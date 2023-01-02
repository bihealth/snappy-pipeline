# -*- coding: utf-8 -*-
"""Wrapper for cnvkit.py access
"""

from snakemake.shell import shell

__author__ = "Manuel Holtgrewe"
__email__ = "manuel.holtgrewe@bih-charite.de"

config = snakemake.config["step_config"][snakemake.config["pipeline_step"]["name"]]["cnvkit"]

exclude = " -- exclude " + " -x ".join(config["exclude"]) if config["exclude"] else ""

shell(
    r"""
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

set -x

# -----------------------------------------------------------------------------

if [[ -n "{config[access]}" ]]
then
    ln -sr {config[access]} {snakemake.output.access}
else
    cnvkit.py access \
        -o {snakemake.output.access} \
        --min-gap-size {config[min_gap_size]} {exclude} \
        {snakemake.config[static_data_config][reference][path]}
fi

fn=$(basename "{snakemake.output.access}")
d=$(dirname "{snakemake.output.access}")
pushd $d
md5sum $fn > $fn.md5
popd
"""
)

# Compute MD5 sums of logs.
shell(
    r"""
md5sum {snakemake.log.log} >{snakemake.log.log_md5}
"""
)

