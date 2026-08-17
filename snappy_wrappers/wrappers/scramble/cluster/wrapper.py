"""CUBI+Snakemake wrapper code for scramble (cluster): Snakemake wrapper.py"""

from typing import TYPE_CHECKING

from snappy_wrappers.snappy_wrapper import ShellWrapper

if TYPE_CHECKING:
    from snakemake.iocontainers import snakemake


ShellWrapper(snakemake).run(
    r"""
# Create out dir
mkdir -p $(dirname {snakemake.output.txt})

# Call tool
cluster_identifier {snakemake.input} > {snakemake.output.txt}
"""
)

