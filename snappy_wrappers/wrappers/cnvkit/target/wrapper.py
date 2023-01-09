# -*- coding: utf-8 -*-
"""Wrapper for cnvkit.py target
"""

from snakemake.shell import shell

__author__ = "Manuel Holtgrewe"
__email__ = "manuel.holtgrewe@bih-charite.de"

step = snakemake.config["pipeline_step"]["name"]
config = snakemake.config["step_config"][step]["cnvkit"]

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
        {snakemake.config[static_data_config][reference][path]}
}}

autobin()
{{
    cnvkit.py autobin --method $1 \
        --fasta {snakemake.config[static_data_config][reference][path]} \
        --access $2 \
        $(if [[ {config[bp_per_bin]} -gt 0 ]]; then \
            echo --bp-per-bin {config[bp_per_bin]}
        fi) \
        $(if [[ {config[target_min_size]} -gt 0 ]]; then \
            echo --target-min-size {config[target_min_size]}
        fi) \
        $(if [[ {config[target_max_size]} -gt 0 ]]; then \
            echo --target-max-size {config[target_max_size]}
        fi) \
        --target-output-bed $tmpdir/target.bed --antitarget-output-bed $tmpdir/antitarget.bed \
        {bams} > $tmpdir/autobin.txt
}}

# -----------------------------------------------------------------------------

target="{config[path_target_regions]}"
target_avg_size={config[target_avg_size]}

if [[ -z "$target" ]] && [[ $target_avg_size -eq 0 ]]
then
    tmpdir=$(mktemp -d $TMPDIR)

    if [[ -n "{bams}" ]]
    then
        if [[ -z "{config[access]}" ]]
        then
            access
            autobin wgs $tmpdir/access.bed

            target=$tmpdir/access.bed
            target_avg_size=$(cat $tmpdir/autobin.txt | grep "Target:" | cut -f 3)
        else
            autobin amplicon "{config[access]}"

            target="{config[access]}"
            target_avg_size=$(cat $tmpdir/autobin.txt | grep "Target:" | cut -f 3)
        fi
    else
        access

        target=$tmpdir/access.bed
        target_avg_size=50000
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
    $(if [[ {config[target_avg_size]} -gt 0 ]]; then \
        echo --avg-size {config[target_avg_size]}
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
