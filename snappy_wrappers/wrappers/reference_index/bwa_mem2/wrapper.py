from typing import TYPE_CHECKING

from snappy_wrappers.snappy_wrapper import ShellWrapper

if TYPE_CHECKING:
    from snakemake.iocontainers import snakemake

extra_args = snakemake.params.args.get("extra_args", "")

ShellWrapper(snakemake).run(
    f"""
bwa-mem2 index {extra_args} -p "$(dirname {{snakemake.output.amb}})/reference" {{snakemake.input.reference}}
"""
)

