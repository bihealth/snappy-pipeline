"""CUBI+Snakemake wrapper code for sequenza (sequenza-utils, pileups)"""

import os
import sys
from typing import TYPE_CHECKING

# The following is required for being able to import snappy_wrappers modules
# inside wrappers.  These run in an "inner" snakemake process which uses its
# own conda environment which cannot see the snappy_pipeline installation.
base_dir = os.path.normpath(os.path.join(os.path.dirname(__file__), "..", "..", "..", ".."))
sys.path.insert(0, base_dir)

from snappy_wrappers.snappy_wrapper import ShellWrapper  # noqa: E402

from snappy_wrappers.tools.genome_windows import yield_contigs  # noqa: E402

__author__ = "Eric Blanc <eric.blanc@bih-charite.de>"

if TYPE_CHECKING:
    from snakemake.iocontainers import snakemake

args = getattr(snakemake.params, "args", {})

genome = args["reference"]
length = args["length"]

f = open(genome + ".fai", "rt")
contigs = " ".join(yield_contigs(f, args.get("ignore_chroms")))
f.close()

extra_arguments = " ".join(
    ["--{} {}".format(k, v) for k, v in args.get("extra_arguments", {}).items()]
)

ShellWrapper(snakemake).run(
    r"""
sequenza-utils bam2seqz \
    -gc {snakemake.input.gc} --fasta {genome} \
    -n {snakemake.input.normal_bam} --tumor {snakemake.input.tumor_bam} \
    -C {contigs} {extra_arguments} \
    | sequenza-utils seqz_binning -s - \
        -w {length} -o {snakemake.output.seqz}
"""
)

