from typing import TYPE_CHECKING

from snappy_wrappers.snappy_wrapper import ShellWrapper

if TYPE_CHECKING:
    from snakemake.iocontainers import snakemake

# samtools faidx writes the .fai next to the FASTA, which is the declared output
ShellWrapper(snakemake).run(
    r"""
samtools faidx {snakemake.input.fasta}
samtools dict -o {snakemake.output.dict} {snakemake.input.fasta}
"""
)
