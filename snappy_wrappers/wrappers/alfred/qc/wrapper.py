from typing import TYPE_CHECKING

from snappy_wrappers.snappy_wrapper import ShellWrapper

if TYPE_CHECKING:
    from snakemake.iocontainers import snakemake

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

reference_path = snakemake.input.reference
genome_path = snakemake.input.reference_genome
targets_bed = snakemake.input.targets_bed

ShellWrapper(snakemake).run(
    r"""
set -x

# Get sorted targets BED file.
zcat --force {targets_bed} \
| awk -F $'\t' 'BEGIN {{ OFS = FS; }} ($2 < $3) {{ print; }}' \
> $TMPDIR/targets.tmp.bed

bedtools sort \
    -i $TMPDIR/targets.tmp.bed \
    -faidx {genome_path} \
| uniq \
> $TMPDIR/targets.bed

# Run "alfred qc".
alfred qc \
    --ignore \
    --reference {reference_path} \
    --bed $TMPDIR/targets.bed \
    --jsonout {snakemake.output.json} \
    --input-file {snakemake.input.bam}
"""
)
