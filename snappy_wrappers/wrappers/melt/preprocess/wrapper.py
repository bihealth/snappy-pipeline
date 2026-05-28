from snappy_wrappers.snappy_wrapper import ShellWrapper

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

args = getattr(snakemake.params, "args", {})
melt_config = args["config"]

ShellWrapper(snakemake).run(
    r"""

ln -sr {snakemake.input.bam} {snakemake.output.orig_bam}
ln -sr {snakemake.input.bai} {snakemake.output.orig_bai}

JAR={melt_config[jar_file]}

java -Xmx13G -jar $JAR \
    Preprocess \
    -bamfile {snakemake.output.orig_bam} \
    -h {snakemake.input.reference}
"""
)
