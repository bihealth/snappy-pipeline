from snappy_wrappers.snappy_wrapper import ShellWrapper
from snappy_wrappers.tools.java import java_heap_mb

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

args = getattr(snakemake.params, "args", {})
melt_arg_exome = "-exome" if args.get("exome", False) else ""
exclude_contigs = "/".join(args["exclude_contigs"])

mem_jvm = java_heap_mb(snakemake)

ShellWrapper(snakemake).run(
    r"""
JAR={snakemake.input.jar}
ME_REFS={snakemake.input.me_refs}
ME_INFIX={args[me_refs_infix]}

mkdir -p {snakemake.output.out}

java -Xmx{mem_jvm}m -jar $JAR \
    IndivAnalysis \
    -b {exclude_contigs} \
    {melt_arg_exome} \
    -h {snakemake.input.reference} \
    -t $ME_REFS/$ME_INFIX/{args[me_type]}_MELT.zip \
    -w {snakemake.output.out} \
    -r {args[read_length]} \
    {extra_arguments} \
    -bamfile {snakemake.input.orig_bam}
"""
)
