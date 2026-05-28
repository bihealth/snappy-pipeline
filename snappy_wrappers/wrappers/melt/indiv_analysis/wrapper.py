from snappy_wrappers.snappy_wrapper import ShellWrapper

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

args = getattr(snakemake.params, "args", {})
melt_config = args["config"]
melt_arg_exome = "-exome" if args.get("exome", False) else ""

ShellWrapper(snakemake).run(
    r"""

JAR={melt_config[jar_file]}
ME_REFS={melt_config[me_refs_path]}
ME_INFIX={melt_config[me_refs_infix]}

java -Xmx13G -jar $JAR \
    IndivAnalysis \
    -b hs37d5/NC_007605 \
    {melt_arg_exome} \
    -h {snakemake.input.reference} \
    -t $ME_REFS/$ME_INFIX/{args[me_type]}_MELT.zip \
    -w $(dirname {snakemake.output.done}) \
    -r 150 \
    -bamfile {snakemake.input.orig_bam}
"""
)
