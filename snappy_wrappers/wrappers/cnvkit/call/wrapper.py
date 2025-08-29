# -*- coding: utf-8 -*-
"""Wrapper vor cnvkit.py call"""

from snakemake.shell import shell

__author__ = "Manuel Holtgrewe"
__email__ = "manuel.holtgrewe@bih-charite.de"

step = snakemake.config["pipeline_step"]["name"]
config = snakemake.config["step_config"][step]["cnvkit"]

center = config["center"]
if center:
    if center in set("mean", "median", "mode", "biweight"):
        center = " --center " + center
    else:
        center = " --center-at" + center
else:
    center = ""

gender = " --gender {}".format(config["gender"]) if config["gender"] else ""
male = " --male-reference" if config["male_reference"] else ""
purity = " --purity {}".foramt(config["purity"]) if config["purity"] > 0 else ""

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

cnvkit.py call \
    --output {snakemake.output.calls} \
    --method {config[calling_method]} \
    --thresholds={config[call_thresholds]} \
    $(if [[ "{config[filter]}" ]]; then \
        echo --filter {config[filter]}
    fi) \
    {center} {gender} {male} \
    --ploidy {config[ploidy]} {purity} \
    {snakemake.input}

d=$(dirname "{snakemake.output.calls}")
pushd $d
fn=$(basename "{snakemake.output.calls}")
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
