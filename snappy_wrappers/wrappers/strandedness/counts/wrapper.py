"""The column of STAR's gene counts (``ReadsPerGene.out.tab``) that the strandedness decision selects.

STAR counts each gene for unstranded, forward and reverse reads in the columns after the gene name.
The decision of the ``rseqc`` wrapper (0, 1 or 2) selects one of them.
"""

import json
from typing import TYPE_CHECKING

from snappy_wrappers.snappy_wrapper import PythonWrapper

if TYPE_CHECKING:
    from snakemake.iocontainers import snakemake

#: The column of the counts of each decision, after the gene name
COLUMNS = {"0": 1, "1": 2, "2": 3}


def main():
    with open(snakemake.input.decision) as f:
        column = COLUMNS[json.load(f)["decision"]]
    with open(snakemake.input.counts, "rb") as src, open(snakemake.output.counts, "wb") as dst:
        for line in src:
            fields = line.rstrip(b"\n").split(b"\t")
            dst.write(fields[0] + b"\t" + fields[column] + b"\n")


PythonWrapper(snakemake).run(main)
