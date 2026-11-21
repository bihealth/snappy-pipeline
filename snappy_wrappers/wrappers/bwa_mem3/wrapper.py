from typing import TYPE_CHECKING

from snappy_wrappers.tools import bwa_mem

if TYPE_CHECKING:
    from snakemake.iocontainers import snakemake

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

bwa_mem.run(snakemake, "bwa-mem3")
