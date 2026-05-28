from typing import TYPE_CHECKING

from snappy_wrappers.snappy_wrapper import ShellWrapper

if TYPE_CHECKING:
    from snakemake.iocontainers import snakemake

extra_args = snakemake.params.args.get("extra_args", "")

ShellWrapper(snakemake).run(
    f"""
minimap2 -t {{snakemake.threads}} -d {{snakemake.output.mmi}} {extra_args} {{snakemake.input.reference}}
"""
)

