# -*- coding: utf-8 -*-
"""Wrapper for cnvkit.py antitarget"""

from snakemake.shell import shell

__author__ = "Manuel Holtgrewe"
__email__ = "manuel.holtgrewe@bih-charite.de"

args = getattr(snakemake.params, "args", {})

target = getattr(snakemake.input, "target", "")
if access := getattr(snakemake.input, "access", ""):
    access = f"--access {access}"

# Avoid compariing with None
args["min_size"] = args["min_size"] if args.get("min_size", None) else 0
args["avg_size"] = args["avg_size"] if args.get("avg_size", None) else 0

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

if [[ -n "{target}" ]]
then
    cnvkit.py antitarget \
        --output {snakemake.output.antitarget} \
        {access} \
        $(if [[ {args[avg_size]} -gt 0 ]]; then \
            echo --avg-size {args[avg_size]}
        fi) \
        $(if [[ {args[min_size]} -gt 0 ]]; then \
            echo --min-size {args[min_size]}
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
