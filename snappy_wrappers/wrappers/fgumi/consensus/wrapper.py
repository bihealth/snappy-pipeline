from typing import TYPE_CHECKING

from snappy_wrappers.tools import fgumi_consensus

if TYPE_CHECKING:
    from snakemake.iocontainers import snakemake

fgumi_consensus.run(snakemake)
