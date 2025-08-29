# -*- coding: utf-8 -*-
"""Wrapper for cnvkit.py plot"""

from snakemake.shell import shell

__author__ = "Manuel Holtgrewe"
__email__ = "manuel.holtgrewe@bih-charite.de"

step = snakemake.config["pipeline_step"]["name"]
config = snakemake.config["step_config"][step]["cnvkit"]

gender = " --gender {}".format(config["gender"]) if config["gender"] else ""
male = " --male-reference" if config["male_reference"] else ""

scatters = [
    snakemake.output.get(x)
    for x in filter(
        lambda x: x.startswith("scatter_chr") and not x.endswith("_md5"), snakemake.output.keys()
    )
]

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

md5()
{{
    d=$(dirname $1)
    f=$(basename $1)
    pushd $d
    md5sum $f > $f.md5
    popd
}}

# -----------------------------------------------------------------------------

unset DISPLAY

if [[ -n "{snakemake.output.diagram}" ]]
then
    cnvkit.py diagram \
        --output {snakemake.output.diagram} \
        --segment {snakemake.input.cns} \
        {gender} {male} \
        --threshold {config[diagram_threshold]} --min-probes {config[diagram_min_probes]} \
        $(if [[ "{config[shift_xy]}" = "False" ]]; then \
            echo --no-shift-xy
        fi) \
        {snakemake.input.cnr}
else
    touch {snakemake.output.diagram}
fi
md5 {snakemake.output.diagram}

if [[ -n "{snakemake.output.scatter}" ]]
then
    cnvkit.py scatter \
        --output {snakemake.output.scatter} \
        --segment {snakemake.input.cns} \
        {gender} {male} \
        {snakemake.input.cnr}
else
    touch {snakemake.output.scatter}
fi
md5 {snakemake.output.scatter}

for scatter in {scatters}
do
    if [[ -n "$scatter" ]]
    then
        chrom=$(basename -s .png $scatter | sed -re "s/.*\.(chr([0-9]+|[XY]))$/\1/")
        cnvkit.py scatter \
            --output $scatter \
            --chromosome $chrom \
            --segment {snakemake.input.cns} \
            {gender} {male} \
            {snakemake.input.cnr}
    else
        touch $scatter
    fi
    md5 $scatter
done
"""
)

# Compute MD5 sums of logs.
shell(
    r"""
md5sum {snakemake.log.log} >{snakemake.log.log_md5}
"""
)
