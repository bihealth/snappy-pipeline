# -*- coding: utf-8 -*-
"""Wrapper for cnvkit.py target"""

from snakemake.shell import shell

__author__ = "Manuel Holtgrewe"
__email__ = "manuel.holtgrewe@bih-charite.de"

args = getattr(snakemake.params, "args", {})
config = args.get("config", {})

bams = " ".join(snakemake.input.get("bams", [""]))

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
        {snakemake.input.reference}
}}

# -----------------------------------------------------------------------------

target="{config[path_target_regions]}"
target_avg_size={config[target_avg_size]}

if [[ -z "$target" ]] && [[ $target_avg_size -eq 0 ]]
then
    tmpdir=$(mktemp -d)

    if [[ -n "{bams}" ]]
    then
        access
        cnvkit.py autobin --method wgs \
            --fasta {snakemake.input.reference} \
            --access $tmpdir/access.bed \
            --bp-per-bin {config[bp_per_bin]} \
            --target-output-bed $tmpdir/target.bed --antitarget-output-bed $tmpdir/antitarget.bed \
            {bams} > $tmpdir/autobin.txt
        target_avg_size=$(cat $tmpdir/autobin.txt | grep "Target:" | cut -f 3)

        if [[ -z "{config[access]}" ]]
        then
            target=$tmpdir/access.bed
        else
            target="{config[access]}"
        fi
    else
        if [[ -z "{config[access]}" ]]
        then
            access
            target=$tmpdir/access.bed
        else
            target="{config[access]}"
        fi
        target_avg_size=5000
    fi
fi

cnvkit.py target \
    --output {snakemake.output.target} \
    $(if [[ -n "{config[annotate]}" ]]; then \
        echo --short-names --annotate {config[annotate]}
    fi) \
    $(if [[ "{config[split]}" = "True" ]]; then \
        echo --split
    fi) \
    $(if [[ $target_avg_size -gt 0 ]]; then \
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
