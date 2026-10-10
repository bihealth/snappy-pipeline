from snappy_wrappers.snappy_wrapper import ShellWrapper

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

args = getattr(snakemake.params, "args", {})

ShellWrapper(snakemake).run(
    r"""

JAR={snakemake.input.jar}
ME_REFS={snakemake.input.me_refs}
ME_INFIX={args[me_refs_infix]}

java -Xmx13G -jar $JAR \
    Genotype \
    -h {snakemake.input.reference} \
    -bamfile {snakemake.input.bam} \
    -p $(dirname {snakemake.input.done}) \
    -t $ME_REFS/$ME_INFIX/{args[me_type]}_MELT.zip \
    -w $(dirname {snakemake.output.done})

"""
)
