from math import ceil

from snappy_wrappers.snappy_wrapper import ShellWrapper

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

args = getattr(snakemake.params, "args", {})
melt_arg_exome = "-exome" if args.get("exome", False) else ""

# As an estimate, give 75% of total memory as usable to JVM
mem_jvm = ceil(snakemake.resources.mem_mb * 0.75)

ShellWrapper(snakemake).run(
    r"""

JAR={snakemake.input.jar}
ME_REFS={snakemake.input.me_refs}
ME_INFIX={args[me_refs_infix]}

java -Xmx{mem_jvm}m -jar $JAR \
    IndivAnalysis \
    -b hs37d5/NC_007605 \
    {melt_arg_exome} \
    -h {snakemake.input.reference} \
    -t $ME_REFS/$ME_INFIX/{args[me_type]}_MELT.zip \
    -w $(dirname {snakemake.output.done}) \
    -r 150 \
    {extra_arguments} \
    -bamfile {snakemake.input.orig_bam}
"""
)
