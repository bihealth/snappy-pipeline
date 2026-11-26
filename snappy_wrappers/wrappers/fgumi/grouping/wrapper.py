from typing import TYPE_CHECKING

from snappy_wrappers.tools import fgumi_grouping

if TYPE_CHECKING:
    from snakemake.iocontainers import snakemake

fgumi_grouping.run(snakemake)
