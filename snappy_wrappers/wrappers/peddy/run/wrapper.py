from typing import TYPE_CHECKING

from snappy_wrappers.snappy_wrapper import ShellWrapper

if TYPE_CHECKING:
    from snakemake.iocontainers import snakemake

ShellWrapper(snakemake).run(
    r"""

python -m peddy \
    --prefix $(dirname {snakemake.output.html})/$(basename {snakemake.output.html} .html) \
    {snakemake.input.vcf} \
    {snakemake.input.ped}
"""
)
