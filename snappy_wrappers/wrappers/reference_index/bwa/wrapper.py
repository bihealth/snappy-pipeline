from typing import TYPE_CHECKING

from snappy_wrappers.snappy_wrapper import ShellWrapper

if TYPE_CHECKING:
    from snakemake.iocontainers import snakemake

algorithm = snakemake.params.args.get("algorithm", "bwtsw")

ShellWrapper(snakemake).run(
    f"""
bwa index -a {algorithm} -p "$(dirname {{snakemake.output.amb}})/reference" {{snakemake.input.reference}}
"""
)

