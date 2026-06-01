from typing import TYPE_CHECKING

from snappy_wrappers.snappy_wrapper import ShellWrapper

if TYPE_CHECKING:
    from snakemake.iocontainers import snakemake

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

args = getattr(snakemake.params, "args", {})
config = args.get("config", {})

ShellWrapper(snakemake).run(
    r"""
# Run actual tools --------------------------------------------------------------------------------

vep --verbose \
    --fasta {snakemake.input.reference} \
    --input_file {snakemake.input.vcf} \
    --output_file {snakemake.output.vcf} \
    --compress_output bgzip \
    --vcf \
    --symbol \
    --terms SO \
    --hgvs \
    --cache \
    --offline \
    --{config[tx_flag]} \
    --force_overwrite \
    --buffer_size {config[buffer_size]} \
    $(if [[ ! -z "{config[cache_dir]}" ]]; then \
        echo --dir_cache {config[cache_dir]}; \
    fi) \
    --cache_version {config[cache_version]} \
    --assembly {config[assembly]} \
    --fork {config[num_threads]} \
    {config[more_flags]}

tabix -f {snakemake.output.vcf}
"""
)
