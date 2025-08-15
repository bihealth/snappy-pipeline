# -*- coding: utf-8 -*-
"""Wrapper vor cnvkit.py segment"""

from snakemake.shell import shell

__author__ = "Manuel Holtgrewe"
__email__ = "manuel.holtgrewe@bih-charite.de"

args = getattr(snakemake.params, "args", {})

method = args["method"]
if method == "cbs" and args["smooth_cbs"]:
    method += " --smooth-cbs"

if float(args["threshold"]) > 0:
    threshold = " --threshold " + str(args["threshold"])
else:
    threshold = ""

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

cnvkit.py segment \
    --output {snakemake.output.segments} \
    --method {method} \
    $(if [[ "{args[drop_low_coverage]}" = "True" ]]; then \
        echo --drop-low-coverage
    fi) \
    {threshold} \
    --drop-outliers {args[drop_outliers]} \
    {snakemake.input}

d=$(dirname "{snakemake.output.segments}")
pushd $d
fn=$(basename "{snakemake.output.segments}")
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
