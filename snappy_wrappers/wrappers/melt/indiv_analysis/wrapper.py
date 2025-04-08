from snakemake.shell import shell

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

melt_config = getattr(snakemake.params, "args", {})
melt_arg_exome = {"sv_calling_targeted": "-exome"}.get(snakemake.params.step_key, "")

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
    IndivAnalysis \
    -b hs37d5/NC_007605 \
    {melt_arg_exome} \
    -h {snakemake.input.reference} \
    -t $ME_REFS/$ME_INFIX/{snakemake.wildcards.me_type}_MELT.zip \
    -w $(dirname {snakemake.output.done}) \
    -r 150 \
    -bamfile {snakemake.input.orig_bam}
"""
)
