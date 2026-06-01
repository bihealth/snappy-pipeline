"""CUBI+Snakemake wrapper code for sequenza GC reference file"""

from typing import TYPE_CHECKING

from snappy_wrappers.snappy_wrapper import ShellWrapper

if TYPE_CHECKING:
    from snakemake.iocontainers import snakemake

__author__ = "Eric Blanc <eric.blanc@bih-charite.de>"

args = getattr(snakemake.params, "args", {})

genome = args["reference"]
length = args["length"]

ShellWrapper(snakemake).run(
    r"""
sequenza-utils gc_wiggle --fasta {genome} -w {length} -o {snakemake.output.gc}
"""
)
