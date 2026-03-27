# -*- coding: utf-8 -*-
"""Wrapper for running VEP variant annotation"""

from snakemake.shell import shell

__author__ = "Eric Blanc"
__email__ = "eric.blanc@bih-charite.de"

args = getattr(snakemake.params, "args", {})
vep_config = args["config"]

# Get shortcuts to step configuration
pick_order = ",".join(vep_config["pick_order"])
script_output_options = " ".join(["--" + x for x in vep_config["output_options"]])
if vep_config["plugins"]:
    plugins = " ".join(["--plugin " + x for x in vep_config["plugins"]])
    if not vep_config["plugins_dir"]:
        raise Exception("Please provide plugins directory if you want to use plugins")
    else:
        plugins_dir = "--dir_plugins " + vep_config["plugins_dir"]
else:
    plugins = ""
    plugins_dir = ""

full = snakemake.output.full if "full" in snakemake.output.keys() else ""

shell(
    r"""
set -x

conda list >{snakemake.log.conda_list}
conda info >{snakemake.log.conda_info}
md5sum {snakemake.log.conda_list} | sed -re "s/  (\.?.+\/)([^\/]+)$/  \2/" > {snakemake.log.conda_list}.md5
md5sum {snakemake.log.conda_info} | sed -re "s/  (\.?.+\/)([^\/]+)$/  \2/" > {snakemake.log.conda_info}.md5

# Also pipe stderr to log file
if [[ -n "{snakemake.log.log}" ]]; then
    if [[ "$(set +e; tty; set -e)" != "" ]]; then
        rm -f "{snakemake.log.log}" && mkdir -p $(dirname {snakemake.log.log})
        exec 2> >(tee -a "{snakemake.log.log}" >&2)
    else
        rm -f "{snakemake.log.log}" && mkdir -p $(dirname {snakemake.log.log})
        echo "No tty, logging disabled" >"{snakemake.log.log}"
    fi
fi

if [[ -n "{full}" ]]
then
    vep --verbose --force_overwrite --offline --cache \
        --fork {vep_config[num_threads]} --buffer_size {vep_config[buffer_size]} \
        --species {vep_config[species]} --cache_version {vep_config[cache_version]} --assembly {vep_config[assembly]} \
        $(if [[ ! -z "{vep_config[cache_dir]}" ]]; then \
            echo --dir_cache {vep_config[cache_dir]}
        fi) \
        {script_output_options} \
        {plugins} \
        {plugins_dir} \
        --{vep_config[tx_flag]} \
        --fasta {snakemake.input.reference} \
        --input_file {snakemake.input.vcf} --format vcf \
        --output_file {full} --vcf --compress_output bgzip
    tabix {full}

    pushd $(dirname {full})
    f=$(basename {full})
    md5sum $f > $f.md5
    md5sum $f.tbi > $f.tbi.md5
    popd
fi

vep --verbose --force_overwrite --offline --cache \
    --fork {vep_config[num_threads]} --buffer_size {vep_config[buffer_size]} \
    --species {vep_config[species]} --cache_version {vep_config[cache_version]} --assembly {vep_config[assembly]} \
    $(if [[ ! -z "{vep_config[cache_dir]}" ]]; then \
        echo --dir_cache {vep_config[cache_dir]}
    fi) \
    {script_output_options} \
    {plugins} \
    {plugins_dir} \
    --pick --pick_order {pick_order} \
    --{vep_config[tx_flag]} \
    --fasta {snakemake.input.reference} \
    --input_file {snakemake.input.vcf} --format vcf \
    --output_file {snakemake.output.vcf} --vcf --compress_output bgzip
tabix {snakemake.output.vcf}

pushd $(dirname {snakemake.output.vcf})
f=$(basename {snakemake.output.vcf})
md5sum $f > $f.md5
md5sum $f.tbi > $f.tbi.md5
popd
"""
)

# Compute MD5 sums of logs.
shell(
    r"""
sleep 1s  # try to wait for log file flush
md5sum {snakemake.log.log} >{snakemake.log.log_md5}
"""
)
