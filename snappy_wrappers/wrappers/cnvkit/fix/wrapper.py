# -*- coding: utf-8 -*-
"""Wrapper for cnvkit.py fix"""

from snakemake.shell import shell

__author__ = "Manuel Holtgrewe"
__email__ = "manuel.holtgrewe@bih-charite.de"

step = snakemake.config["pipeline_step"]["name"]
config = snakemake.config["step_config"][step]["cnvkit"]

if "ref" in snakemake.input.keys():
    ref = snakemake.input.target
elif "path_panel_of_normals" in config.keys():
    ref = config["path_panel_of_normals"]
else:
    raise Exception("Unsupported naming")

gender = " --gender {}".format(config["gender"]) if config["gender"] else ""
male = " --male-reference" if config["male_reference"] else ""
no_gc = " --no-gc" if not config["gc_correction"] else ""
no_edge = " --no-edge" if not config["edge_correction"] else ""
no_rmask = " --no-rmask" if not config["rmask_correction"] else ""

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

cnvkit.py fix \
    --output {snakemake.output.ratios} \
    {gender} {male} {no_gc} {no_edge} {no_rmask} \
    {snakemake.input.target} \
    {snakemake.input.antitarget} \
    {ref}

d=$(dirname "{snakemake.output.ratios}")
pushd $d
fn=$(basename "{snakemake.output.ratios}")
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
