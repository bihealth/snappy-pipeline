from math import ceil

from snappy_wrappers.snappy_wrapper import ShellWrapper

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

args = getattr(snakemake.params, "args", {})

# As an estimate, give 75% of total memory as usable to JVM
mem_jvm = ceil(snakemake.resources.mem_mb * 0.75)

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
    -w $(dirname {snakemake.output.done})

"""
)
