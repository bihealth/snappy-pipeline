# -*- coding: utf-8 -*-
"""Wrapper for cnvkit.py report
"""

from snakemake.shell import shell

__author__ = "Manuel Holtgrewe"
__email__ = "manuel.holtgrewe@bih-charite.de"

step = snakemake.config["pipeline_step"]["name"]
config = snakemake.config["step_config"][step]["cnvkit"]

for param in (
    "breaks_min_probes",
    "genemetrics_min_probes",
    "genemetrics_threshold",
    "drop_low_coverage",
    "genemetrics_alpha",
    "genemetrics_bootstrap",
    "segmetrics_alpha",
    "segmetrics_bootstrap",
    "smooth_bootstrap",
):
    if not param in config.keys():
        config[param] = ""

gender = " --gender {}".format(config["gender"]) if config["gender"] else ""
male = " --male-reference" if config["male_reference"] else ""

input_target = snakemake.input.get("target", "")
input_antitarget = snakemake.input.get("antitarget", "")
input_cnr = snakemake.input.get("cnr", "")
input_cns = snakemake.input.get("cns", "")

output_breaks = snakemake.output.get("breaks", "")
output_genemetrics = snakemake.output.get("genemetrics", "")

output_segmetrics = snakemake.output.get("segmetrics", "")
output_sex = snakemake.output.get("sex", "")
output_metrics = snakemake.output.get("metrics", "")

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

if [[ -n "{output_breaks}" ]]
then
    if [[ -n "{input_cnr}" ]] && [[ -n "{input_cns}" ]]
    then
        cnvkit.py breaks \
            --output {output_breaks} \
            --min-probes {config[breaks_min_probes]} \
            {input_cnr} {input_cns}
    else
        touch {output_breaks}
    fi
    md5 {output_breaks}
fi

if [[ -n "{output_genemetrics}" ]]
then
    if [[ -n "{input_cnr}" ]] && [[ -n "{input_cns}" ]]
    then
        cnvkit.py genemetrics \
            --output {output_genemetrics} \
            --segment {input_cns} \
            --min-probes {config[genemetrics_min_probes]} --threshold {config[genemetrics_threshold]} \
            $(if [[ "{config[drop_low_coverage]}" = "True" ]]; then \
                echo --drop-low-coverage
            fi) \
            {gender} {male} \
            --mean --median --mode --ttest --stdev --sem --mad --mse --iqr --bivar --ci --pi \
            --alpha {config[genemetrics_alpha]} --bootstrap {config[genemetrics_bootstrap]} \
            {input_cnr}
    else
        touch {output_genemetrics}
    fi
    md5 {output_genemetrics}
fi

if [[ -n "{output_segmetrics}" ]]
then
    if [[ -n "{input_cnr}" ]] && [[ -n "{input_cns}" ]]
    then
        cnvkit.py segmetrics \
            --output {output_segmetrics} \
            --segment {input_cns} \
            --mean --median --mode --t-test --stdev --sem --mad --mse --iqr --bivar --ci --pi \
            --alpha {config[segmetrics_alpha]} --bootstrap {config[segmetrics_bootstrap]} \
            $(if [[ "{config[smooth_bootstrap]}" = "True" ]]; then \
                echo --smooth-bootstrap
            fi) \
            $(if [[ "{config[drop_low_coverage]}" = "True" ]]; then \
                echo --drop-low-coverage
            fi) \
            {input_cnr}
    else
        touch {output_segmetrics}
    fi
    md5 {output_segmetrics}
fi

if [[ -n "{output_sex}" ]]
then
    if [[ -n "{input_target}" || -n "{input_antitarget}" || -n "{input_cnr}" || -n "{input_cns}" ]]
    then
        cnvkit.py sex \
            --output {output_sex} \
            $(if [[ "{male}" = "True" ]]; then \
                echo --male-reference
            fi) \
            $(if [[ -n "{input_target}" ]]; then \
                echo {input_target}
            fi) \
            $(if [[ -n "{input_antitarget}" ]]; then \
                echo {input_antitarget}
            fi) \
            $(if [[ -n "{input_cnr}" ]]; then \
                echo {input_cnr}
            fi) \
            $(if [[ -n "{input_cns}" ]]; then \
                echo {input_cns}
            fi)
    else
        touch {output_sex}
    fi
    md5 {output_sex}
fi

if [[ -n "{output_metrics}" ]]
then
    if [[ -n "{input_target}" || -n "{input_antitarget}" || -n "{input_cnr}" || -n "{input_cns}" ]]
    then
        cnvkit.py metrics \
            --output {output_metrics} \
            $(if [[ "config[drop_low_coverage]" = "True" ]]; then \
                echo --drop-low-coverage
            fi) \
            $(if [[ -n "{input_target}" ]]; then \
                echo {input_target}
            fi) \
            $(if [[ -n "{input_antitarget}" ]]; then \
                echo {input_antitarget}
            fi) \
            $(if [[ -n "{input_cnr}" ]]; then \
                echo {input_cnr}
            fi) \
            $(if [[ -n "{input_cns}" ]]; then \
                echo --segments {input_cns}
            fi)
    else
        touch {output_metrics}
    fi
    md5 {output_metrics}
fi
"""
)

# Compute MD5 sums of logs.
shell(
    r"""
md5sum {snakemake.log.log} >{snakemake.log.log_md5}
"""
)
