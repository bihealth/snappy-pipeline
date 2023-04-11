# -*- coding: utf-8 -*-
"""Wrapper for running VEP variant annotation
"""

from snakemake.shell import shell

__author__ = "Eric Blanc"
__email__ = "eric.blanc@bih-charite.de"

# Get shortcuts to step configuration
current_step = snakemake.config["pipeline_step"]["name"]
vep_config = snakemake.config["step_config"][current_step]["vep"]
script_output_options = " ".join(["--" + x for x in vep_config["output_options"]])

shell(
    r"""
set -x

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

vep --verbose --force_overwrite --offline --cache \
    --fork {vep_config[num_threads]} --buffer_size {vep_config[buffer_size]} \
    --species {vep_config[species]} --cache_version {vep_config[cache_version]} --assembly {vep_config[assembly]} \
    $(if [[ ! -z "{vep_config[cache_dir]}" ]]; then \
        echo --dir_cache {vep_config[cache_dir]}
    fi) \
    {script_output_options} \
    $(if [[ "{vep_config[pick]}" = "yes" ]]; then \
        echo "--pick"
    fi) \
    --{vep_config[tx_flag]} \
    --fasta {snakemake.config[static_data_config][reference][path]} \
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
