"""Duplicate marking with the molecular barcodes of AGeNT (LocatIt), as mbcs ran it"""

from typing import TYPE_CHECKING

from snappy_wrappers.snappy_wrapper import ShellWrapper

if TYPE_CHECKING:
    from snakemake.iocontainers import snakemake

args = snakemake.params.args
input_filter_args = " ".join(args["input_filter_args"])
consensus_filter_args = " ".join(args["consensus_filter_args"])

ShellWrapper(snakemake).run(
    r"""
java {args[java_options]} -Djava.io.tmpdir=$TMPDIR -jar {snakemake.input.jar} {extra_arguments} \
    --bed-file={snakemake.input.baits} \
    -c={args[consensus_mode]} \
    -o={snakemake.output.bam} \
    -f {input_filter_args} -F {consensus_filter_args} \
    {snakemake.input.bam}
samtools index {snakemake.output.bam}
"""
)
