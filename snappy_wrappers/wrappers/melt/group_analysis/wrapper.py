from snakemake.shell import shell

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

args = getattr(snakemake.params, "args", {})
melt_config = args["config"]

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

java -jar -Xmx13G -jar $JAR \
    GroupAnalysis \
    -h {snakemake.input.reference} \
    -t $ME_REFS/$ME_INFIX/{args[me_type]}_MELT.zip \
    $(if [[ $ME_REFS == *37* ]] || [[ $ME_REFS == *hg19* ]]; then
        echo -v $ME_REFS/../../prior_files/{args[me_type]}.1KGP.sites.vcf;
    fi) \
    -w $(dirname {snakemake.output.done}) \
    -r 150 \
    -n {melt_config[genes_file]} \
    -discoverydir $(dirname $(echo {snakemake.input} | cut -d ' ' -f 1))
"""
)
