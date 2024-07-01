# -*- coding: utf-8 -*-
"""Wrapper for cnvkit.py antitarget"""

from snakemake.shell import shell

__author__ = "Manuel Holtgrewe"
__email__ = "manuel.holtgrewe@bih-charite.de"

step = snakemake.config["pipeline_step"]["name"]
config = snakemake.config["step_config"][step]["cnvkit"]

target = snakemake.input.get("target", "")

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

if [[ -n "{config[path_target_regions]}" ]]
then
    cnvkit.py antitarget \
        --output {snakemake.output.antitarget} \
        $(if [[ -n "{config[access]}" ]]; then \
            echo --access {config[access]}
        fi) \
        $(if [[ {config[antitarget_avg_size]} -gt 0 ]]; then \
            echo --avg-size {config[antitarget_avg_size]}
        fi) \
        $(if [[ {config[min_size]} -gt 0 ]]; then \
            echo --min-size {config[min_size]}
        fi) \
        {target}
else
    touch {snakemake.output.antitarget}
fi

fn=$(basename "{snakemake.output.antitarget}")
d=$(dirname "{snakemake.output.antitarget}")
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
