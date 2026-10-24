"""Base quality score recalibration with GATK BaseRecalibrator and ApplyBQSR"""

from typing import TYPE_CHECKING

from snappy_wrappers.snappy_wrapper import ShellWrapper

if TYPE_CHECKING:
    from snakemake.iocontainers import snakemake

ShellWrapper(snakemake).run(
    r"""
gatk BaseRecalibrator \
    --tmp-dir $TMPDIR \
    --input {snakemake.input.bam} \
    --reference {snakemake.input.reference} \
    --known-sites {snakemake.input.common_variants} \
    --output {snakemake.output.recal_table} \
    {extra_arguments}

# samtools writes the index as <bam>.bai, GATK as <name>.bai
gatk ApplyBQSR \
    --tmp-dir $TMPDIR \
    --input {snakemake.input.bam} \
    --reference {snakemake.input.reference} \
    --bqsr-recal-file {snakemake.output.recal_table} \
    --create-output-bam-index false \
    --output {snakemake.output.bam}
samtools index {snakemake.output.bam}
"""
)
