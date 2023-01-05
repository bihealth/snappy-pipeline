from snakemake.shell import shell

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

shell(
    r"""
# -----------------------------------------------------------------------------
# Redirect stderr to log file by default and enable printing executed commands
exec 2> >(tee -a "{snakemake.log.log}")
set -x
# -----------------------------------------------------------------------------

ln -sr {snakemake.input.bam} {snakemake.output.orig_bam}
ln -sr {snakemake.input.bai} {snakemake.output.orig_bai}

JAR={snakemake.config[step_config][sv_calling_targeted][melt][jar_file]}

java -Xmx4G -jar $JAR \
    Preprocess \
    -bamfile {snakemake.output.orig_bam} \
    -h {snakemake.config[static_data_config][reference][path]}
"""
)
