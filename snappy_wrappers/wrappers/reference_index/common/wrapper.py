from typing import TYPE_CHECKING

from snappy_wrappers.snappy_wrapper import ShellWrapper

if TYPE_CHECKING:
    from snakemake.iocontainers import snakemake

ShellWrapper(snakemake).run(
    r"""
samtools faidx {snakemake.input.reference}
cp -f {snakemake.input.reference}.fai {snakemake.output.reference_fai}
samtools dict -o {snakemake.output.reference_dict} {snakemake.input.reference}
cut -f1,2 {snakemake.output.reference_fai} > {snakemake.output.reference_genome}
"""
)

