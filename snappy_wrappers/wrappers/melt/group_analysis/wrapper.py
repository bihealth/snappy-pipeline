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

JAR={melt_config[jar_file]}
ME_REFS={melt_config[me_refs_path]}
ME_INFIX={melt_config[me_refs_infix]}

java -jar -Xmx8G -jar $JAR \
    GroupAnalysis \
    -h {snakemake.config[static_data_config][reference][path]} \
    -t $ME_REFS/$ME_INFIX/{snakemake.wildcards.me_type}_MELT.zip \
    $(if [[ $ME_REFS == *37* ]] || [[ $ME_REFS == *hg19* ]]; then
        echo -v $ME_REFS/../../prior_files/{snakemake.wildcards.me_type}.1KGP.sites.vcf;
    fi) \
    -w $(dirname {snakemake.output.done}) \
    -r 150 \
    -n {melt_config[genes_file]} \
    -discoverydir $(dirname $(echo {snakemake.input} | cut -d ' ' -f 1))
"""
)
