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

# peddy names its pedigree {{prefix}}.peddy.ped, the declared name is {{prefix}}.ped
mv $(dirname {snakemake.output.html})/$(basename {snakemake.output.html} .html).peddy.ped {snakemake.output.ped}
"""
)
