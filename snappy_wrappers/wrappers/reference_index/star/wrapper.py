from typing import TYPE_CHECKING

from snappy_wrappers.snappy_wrapper import ShellWrapper

if TYPE_CHECKING:
    from snakemake.iocontainers import snakemake

features = snakemake.input.get("features", "")

feature_arg = f"--sjdbGTFfile {features}" if features else ""

ShellWrapper(snakemake).run(
    f"""
STAR \\
  --runMode genomeGenerate \\
  --runThreadN {{snakemake.threads}} \\
  --genomeDir "$(dirname {{snakemake.output.index_Genome}})" \\
  --genomeFastaFiles {{snakemake.input.reference}} \\
  {feature_arg} \\
  {{extra_arguments}}
"""
)
