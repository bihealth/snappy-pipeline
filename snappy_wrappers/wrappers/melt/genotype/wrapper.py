from snakemake.shell import shell

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

melt_config = getattr(snakemake.params, "args", {})

shell(
    r"""
# -----------------------------------------------------------------------------
# Redirect stderr to log file by default and enable printing executed commands
exec 2> >(tee -a "{snakemake.log.log}")
set -x
# -----------------------------------------------------------------------------

JAR={melt_config[jar_file]}
ME_REFS={melt_config[me_refs_path]}
ME_INFIX={melt_config[me_refs_infix]}

java -Xmx13G -jar $JAR \
    Genotype \
    -h {melt_config[reference]} \
    -bamfile {snakemake.input.bam} \
    -p $(dirname {snakemake.input.done}) \
    -t $ME_REFS/$ME_INFIX/{snakemake.wildcards.me_type}_MELT.zip \
    -w $(dirname {snakemake.output.done})

"""
)
