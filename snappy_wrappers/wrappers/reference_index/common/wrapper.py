from typing import TYPE_CHECKING

from snappy_wrappers.snappy_wrapper import ShellWrapper

if TYPE_CHECKING:
    from snakemake.iocontainers import snakemake

ShellWrapper(snakemake).run(
    r"""
# Write the indexes here, not next to the input FASTA, which belongs to the upstream task. The
# .gzi of a bgzipped FASTA goes to the temporary directory; no output needs it.
samtools faidx \
    --fai-idx {snakemake.output.reference_fai} \
    --gzi-idx $TMPDIR/reference.gzi \
    {snakemake.input.reference}
samtools dict -o {snakemake.output.reference_dict} {snakemake.input.reference}
cut -f1,2 {snakemake.output.reference_fai} > {snakemake.output.reference_genome}
"""
)
