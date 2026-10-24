from math import ceil

from snappy_wrappers.snappy_wrapper import ShellWrapper

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

# As an estimate, give 75% of total memory as usable to JVM
mem_jvm = ceil(snakemake.resources.mem_mb * 0.75)

ShellWrapper(snakemake).run(
    r"""
ln -sr {snakemake.input.bam} {snakemake.output.orig_bam}
ln -sr {snakemake.input.bai} {snakemake.output.orig_bai}

JAR={snakemake.input.jar}

java -Xmx{mem_jvm}m -jar $JAR \
    Preprocess \
    -bamfile {snakemake.output.orig_bam} \
    -h {snakemake.input.reference}
"""
)
