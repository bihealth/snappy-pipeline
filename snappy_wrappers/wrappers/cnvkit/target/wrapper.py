# -*- coding: utf-8 -*-
"""Wrapper for cnvkit.py target
"""

from snakemake.shell import shell

__author__ = "Manuel Holtgrewe"
__email__ = "manuel.holtgrewe@bih-charite.de"

step = snakemake.config["pipeline_step"]["name"]
config = snakemake.config["step_config"][step]["cnvkit"]

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

cnvkit.py target \
    --output {snakemake.output.target} \
    $(if [[ -n "{config[annotate]}" ]]; then \
        echo --short-names --annotate {config[annotate]}
    fi) \
    $(if [[ "{config[split]}" = "True" ]]; then \
        echo --split
    fi) \
    $(if [[ {config[target_avg_size]} -gt 0 ]]; then \
        echo --avg-size {config[target_avg_size]}
    fi) \
    {config[path_target_regions]}

fn=$(basename "{snakemake.output.target}")
d=$(dirname "{snakemake.output.target}")
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
