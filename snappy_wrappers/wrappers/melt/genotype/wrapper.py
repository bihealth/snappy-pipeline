from snappy_wrappers.snappy_wrapper import ShellWrapper
from snappy_wrappers.tools.java import java_heap_mb

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

args = getattr(snakemake.params, "args", {})

mem_jvm = java_heap_mb(snakemake)

ShellWrapper(snakemake).run(
    r"""
JAR={snakemake.input.jar}
ME_REFS={snakemake.input.me_refs}
ME_INFIX={args[me_refs_infix]}

java -Xmx{mem_jvm}m -jar $JAR \
    Genotype \
    -h {snakemake.input.reference} \
    -bamfile {snakemake.input.bam} \
    -p $(dirname {snakemake.input.pre_geno}) \
    -t $ME_REFS/$ME_INFIX/{args[me_type]}_MELT.zip \
    -w $(dirname {snakemake.output.tsv})
"""
)
