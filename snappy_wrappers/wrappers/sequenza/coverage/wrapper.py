"""CUBI+Snakemake wrapper code for sequenza (sequenza-utils, pileups)"""

from typing import TYPE_CHECKING

from snappy_wrappers.snappy_wrapper import ShellWrapper

from snappy_wrappers.tools.genome_windows import yield_contigs

__author__ = "Eric Blanc <eric.blanc@bih-charite.de>"

if TYPE_CHECKING:
    from snakemake.iocontainers import snakemake

args = getattr(snakemake.params, "args", {})

genome = snakemake.input.reference
length = args["length"]

f = open(snakemake.input.reference_fai, "rt")
contigs = " ".join(yield_contigs(f, args.get("ignore_contigs")))
f.close()

bam2seqz_options = " ".join(
    ["--{} {}".format(k, v) for k, v in args.get("bam2seqz_options", {}).items()]
)

ShellWrapper(snakemake).run(
    r"""
sequenza-utils bam2seqz \
    -gc {snakemake.input.gc} --fasta {genome} \
    -n {snakemake.input.normal_bam} --tumor {snakemake.input.tumor_bam} \
    -C {contigs} {bam2seqz_options} {extra_arguments} \
    | sequenza-utils seqz_binning -s - \
        -w {length} -o {snakemake.output.seqz}
"""
)
