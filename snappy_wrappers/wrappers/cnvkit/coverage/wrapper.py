# -*- coding: utf-8 -*-
"""Wrapper for cnvkit.py coverage"""

from snakemake.shell import shell

__author__ = "Manuel Holtgrewe"
__email__ = "manuel.holtgrewe@bih-charite.de"

step = snakemake.config["pipeline_step"]["name"]
config = snakemake.config["step_config"][step]["cnvkit"]

# During panel_of_normals step, the target regions are created by the target substep.
# During somatic CNV calling (both exome & wgs), the target regions are obtained from the configuration
if "target" in snakemake.input.keys():
    target = snakemake.input.target
elif "path_target" in config.keys():
    target = config["path_target"]
else:
    raise Exception("Unsupported naming")

# Same for antitarget regions
if "antitarget" in snakemake.input.keys():
    antitarget = snakemake.input.antitarget
elif "path_antitarget" in config.keys():
    antitarget = config["path_antitarget"]
else:
    raise Exception("Unsupported naming")

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

# Function definitions ---------------------------------------------------------

coverage()
{{
    cnvkit.py coverage \
        --fasta {snakemake.config[static_data_config][reference][path]} \
        --min-mapq {config[min_mapq]} \
        --processes {snakemake.threads} \
        {snakemake.input.bam} \
        --output $2 $1
}}

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

coverage {target} {snakemake.output.target}
md5 {snakemake.output.target}

coverage {antitarget} {snakemake.output.antitarget}
md5 {snakemake.output.antitarget}
"""
)

# Compute MD5 sums of logs.
shell(
    r"""
md5sum {snakemake.log.log} >{snakemake.log.log_md5}
"""
)
