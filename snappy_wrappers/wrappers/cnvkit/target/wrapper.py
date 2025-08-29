# -*- coding: utf-8 -*-
"""Wrapper for cnvkit.py target"""

from snakemake.shell import shell

__author__ = "Manuel Holtgrewe"
__email__ = "manuel.holtgrewe@bih-charite.de"

args = getattr(snakemake.params, "args", {})

bams = " ".join(getattr(snakemake.input, "bams", [""]))
reference = getattr(snakemake.input, "reference", "")

target = getattr(snakemake.input, "target", "")
access = getattr(snakemake.input, "access", "")

if annotate := getattr(snakemake.input, "annotate", ""):
    annotate = f"--short-names --annotate {annotate}"

# Avoid testing floats in bash
target_avg_size = args["avg_size"] if args.get("avg_size", None) else ""

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

access()
{{
    cnvkit.py access \
        -o $tmpdir/access.bed \
        {reference}
}}

# -----------------------------------------------------------------------------

target="{target}"
target_avg_size={target_avg_size}

# No target file -> in WGS mode -> targets are accessible regions
if [[ -z "$target" ]] && [[ -z "$target_avg_size" ]]
then
    tmpdir=$(mktemp -d)

    # Normals are available -> compute target_avg_size from bams, without using user-defined access file
    if [[ -n "{bams}" ]]
    then
        access
        cnvkit.py autobin --method wgs \
            --fasta {reference} \
            --access $tmpdir/access.bed \
            --bp-per-bin {args[bp_per_bin]} \
            --target-output-bed $tmpdir/target.bed --antitarget-output-bed $tmpdir/antitarget.bed \
            {bams} > $tmpdir/autobin.txt
        target_avg_size=$(cat $tmpdir/autobin.txt | grep "Target:" | cut -f 3)

        if [[ -z "{access}" ]]
        then
            target=$tmpdir/access.bed
        else
            target="{access}"
        fi
    else
        if [[ -z "{access}" ]]
        then
            access
            target=$tmpdir/access.bed
        else
            target="{access}"
        fi
        target_avg_size=5000
    fi
fi

cnvkit.py target \
    --output {snakemake.output.target} \
    {annotate} \
    $(if [[ "{args[split]}" = "True" ]]; then \
        echo --split
    fi) \
    $(if [[ -n "$target_avg_size" ]]; then \
        echo --avg-size $target_avg_size
    fi) \
    $target

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
