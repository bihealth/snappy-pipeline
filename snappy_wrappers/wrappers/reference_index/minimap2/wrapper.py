from typing import TYPE_CHECKING

from snappy_wrappers.snappy_wrapper import ShellWrapper

if TYPE_CHECKING:
    from snakemake.iocontainers import snakemake

ShellWrapper(snakemake).run(
    r"""
minimap2 -t {snakemake.threads} -d {snakemake.output.mmi} {extra_arguments} {snakemake.input.reference}
"""
)

