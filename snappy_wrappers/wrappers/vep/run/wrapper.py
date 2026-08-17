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

dir_cache = f"--dir_cache {vep_config['cache_dir']}" if vep_config.get("cache_dir") else ""
num_threads = getattr(snakemake, "threads", vep_config.get("num_threads", 1))
full = snakemake.output.full if "full" in snakemake.output.keys() else ""

ShellWrapper(snakemake).run(
    r"""
# Helper function to normalize, annotate with VEP, and index output
run_vep() {
    local outfile="$1"
    shift

    bcftools norm --multiallelics -any {snakemake.input.vcf} --threads {num_threads} --force | \
    vep --verbose --force_overwrite --offline --cache \
        --fork {vep_config[num_threads]} --buffer_size {vep_config[buffer_size]} \
        --species {vep_config[species]} --cache_version {vep_config[cache_version]} --assembly {vep_config[assembly]} \
        {dir_cache} \
        {script_output_options} \
        --{vep_config[tx_flag]} \
        --fasta {snakemake.input.reference} \
        --format vcf --vcf --compress_output bgzip \
        --output_file "$outfile" "$@"

    tabix "$outfile"
}

# full annotation if requested
if [[ -n "{full}" ]]; then
    run_vep "{full}"
fi

run_vep "{snakemake.output.vcf}" --pick --pick_order {pick_order}
"""
)
