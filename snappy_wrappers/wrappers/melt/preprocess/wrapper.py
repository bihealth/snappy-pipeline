from snakemake.shell import shell

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

melt_config = snakemake.config["step_config"][snakemake.params.step_key]["melt"]

shell(
    r"""
# -----------------------------------------------------------------------------
# Redirect stderr to log file by default and enable printing executed commands
exec 2> >(tee -a "{snakemake.log.log}")
set -x
# -----------------------------------------------------------------------------

ln -sr {snakemake.input.bam} {snakemake.output.orig_bam}
ln -sr {snakemake.input.bai} {snakemake.output.orig_bai}

JAR={melt_config[jar_file]}

java -Xmx8G -jar $JAR \
    Preprocess \
    -bamfile {snakemake.output.orig_bam} \
    -h {snakemake.config[static_data_config][reference][path]}
"""
)
