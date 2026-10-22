from typing import TYPE_CHECKING

from snappy_wrappers.snappy_wrapper import ShellWrapper

if TYPE_CHECKING:
    from snakemake.iocontainers import snakemake

ShellWrapper(snakemake).run(
    r"""
bwa-mem2 index {extra_arguments} -p "$(dirname {snakemake.output.amb})/reference" {snakemake.input.reference}
"""
)

