# -*- coding: utf-8 -*-
"""Wrapper for cnvkit.py fix"""

from snakemake.shell import shell

__author__ = "Manuel Holtgrewe"
__email__ = "manuel.holtgrewe@bih-charite.de"

args = getattr(snakemake.params, "args", {})

gender = " --gender {}".format(args["gender"]) if args.get("gender", None) else ""
male = " --male-reference" if args.get("male_reference", False) else ""
no_gc = " --no-gc" if not args["gc_correction"] else ""
no_edge = " --no-edge" if not args["edge_correction"] else ""
no_rmask = " --no-rmask" if not args["rmask_correction"] else ""

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
    {snakemake.input.ref}

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
