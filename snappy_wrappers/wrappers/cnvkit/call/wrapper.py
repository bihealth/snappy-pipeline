# -*- coding: utf-8 -*-
"""Wrapper vor cnvkit.py call"""

from snakemake.shell import shell

__author__ = "Manuel Holtgrewe"
__email__ = "manuel.holtgrewe@bih-charite.de"

args = getattr(snakemake.params, "args", {})

if center := args["center"]:
    if center in ("mean", "median", "mode", "biweight"):
        center = " --center " + center
    else:
        center = " --center-at" + str(center)
else:
    center = ""
    
gender = " --gender {}".format(args["gender"]) if args.get("gender", None) else ""
male = " --male-reference" if args.get("male_reference", False) else ""
purity = " --purity {}".format(args["purity"]) if args.get("purity", 0.0) > 0.0 else ""

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
    --method {args[method]} \
    --thresholds={args[thresholds]} \
    $(if [[ -n "{args[filter]}" ]]; then \
        echo --filter {args[filter]}
    fi) \
    {center} {gender} {male} \
    --ploidy {args[ploidy]} {purity} \
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
