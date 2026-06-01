# -*- coding: utf-8 -*-
"""Wrapper for running VEP variant annotation"""

from typing import TYPE_CHECKING

from snappy_wrappers.snappy_wrapper import ShellWrapper

if TYPE_CHECKING:
    from snakemake.iocontainers import snakemake

__author__ = "Eric Blanc"
__email__ = "eric.blanc@bih-charite.de"

args = getattr(snakemake.params, "args", {})
vep_config = args["config"]

# Get shortcuts to step configuration
pick_order = ",".join(vep_config["pick_order"])
script_output_options = " ".join(["--" + x for x in vep_config["output_options"]])

full = snakemake.output.full if "full" in snakemake.output.keys() else ""

ShellWrapper(snakemake).run(
    r"""
if [[ -n "{full}" ]]
then
    vep --verbose --force_overwrite --offline --cache \
        --fork {vep_config[num_threads]} --buffer_size {vep_config[buffer_size]} \
        --species {vep_config[species]} --cache_version {vep_config[cache_version]} --assembly {vep_config[assembly]} \
        $(if [[ ! -z "{vep_config[cache_dir]}" ]]; then \
            echo --dir_cache {vep_config[cache_dir]}
        fi) \
        {script_output_options} \
        --{vep_config[tx_flag]} \
        --fasta {snakemake.input.reference} \
        --input_file {snakemake.input.vcf} --format vcf \
        --output_file {full} --vcf --compress_output bgzip
    tabix {full}
fi

vep --verbose --force_overwrite --offline --cache \
    --fork {vep_config[num_threads]} --buffer_size {vep_config[buffer_size]} \
    --species {vep_config[species]} --cache_version {vep_config[cache_version]} --assembly {vep_config[assembly]} \
    $(if [[ ! -z "{vep_config[cache_dir]}" ]]; then \
        echo --dir_cache {vep_config[cache_dir]}
    fi) \
    {script_output_options} \
    --pick --pick_order {pick_order} \
    --{vep_config[tx_flag]} \
    --fasta {snakemake.input.reference} \
    --input_file {snakemake.input.vcf} --format vcf \
    --output_file {snakemake.output.vcf} --vcf --compress_output bgzip
tabix {snakemake.output.vcf}
"""
)
