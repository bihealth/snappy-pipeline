from typing import TYPE_CHECKING

from snappy_wrappers.snappy_wrapper import ShellWrapper

if TYPE_CHECKING:
    from snakemake.iocontainers import snakemake

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

num_input_files = len(snakemake.input.vcf)

ShellWrapper(snakemake).run(
    r"""
# Run actual tools --------------------------------------------------------------------------------

if [[ {num_input_files} -gt 1 ]]; then
    python $(dirname {__file__})/../../../tools/gcnv_merge_vcfs.py \
        {snakemake.output.vcf} \
        {snakemake.input.vcf}
    tabix -f {snakemake.output.vcf}
else
    cp -a {snakemake.input.vcf} {snakemake.output.vcf}
    cp -a {snakemake.input.vcf}.tbi {snakemake.output.vcf}.tbi
fi
"""
)
