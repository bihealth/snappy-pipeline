# -*- coding: utf-8 -*-
"""Wrapper vor cnvkit.py genome2access
"""

from snakemake.shell import shell

__author__ = "Manuel Holtgrewe"
__email__ = "manuel.holtgrewe@bih-charite.de"

config = snakemake.config["step_config"][snakemake.config["pipeline_step"]["name"]]["cnvkit"]

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

# Function definitions ---------------------------------------------------------
md5()
{{
    set -x

    fn=$1
    f=$(basename $fn)
    d=$(dirname $fn)
    pushd $d
    md5sum $f > $f.md5
    popd
}}

# -----------------------------------------------------------------------------

cnvkit.py autobin --method "hybrid" \
    --fasta {snakemake.config[static_data_config][reference][path]} \
    --access {snakemake.input.access} \
    --targets {config[path_target_regions]} \
    --bp-per-bin {config[bp_per_bin]} \
    --target-min-size {config[target_min_size]} --target-max-size {config[target_max_size]} \
    --antitarget-min-size {config[antitarget_min_size]} --antitarget-max-size {config[antitarget_max_size]} \
    $(if [[ -n "{config[annotate]}" ]] ; then \
        echo --annotate {config[annotate]} --short-names
    fi) \
    --target-output-bed {snakemake.output.target} --antitarget-output-bed {snakemake.output.antitarget} \
    {snakemake.input.bams}

md5 {snakemake.output.target}
md5 {snakemake.output.antitarget}
"""
)
