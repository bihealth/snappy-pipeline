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

JAR={snakemake.config[step_config][sv_calling_targeted][melt][jar_file]}
ME_REFS={melt_config[me_refs_path]}
ME_INFIX={melt_config[me_refs_infix]}

genotype_dir=$(dirname {snakemake.input.genotype} | head -n 1)
ls $genotype_dir/*.{snakemake.wildcards.me_type}.tsv \
| sort \
> {snakemake.output.list_txt}

# TODO: allowing 100% no-call wise?
java -Xmx6G -jar $JAR \
    MakeVCF \
    -genotypingdir $genotype_dir \
    -h {snakemake.config[static_data_config][reference][path]} \
    -j 100 \
    -t $ME_REFS/$ME_INFIX/{snakemake.wildcards.me_type}_MELT.zip \
    -p $(dirname {snakemake.input.group_analysis}) \
    -w $(dirname {snakemake.output.done}) \
    -o $(dirname {snakemake.output.done})

## Make file more VCF conforming XXX TODO XXX
#perl -p -e 's/ID=GL,Number=3/ID=GL,Number=G/' $(dirname {snakemake.output.done})/{snakemake.wildcards.me_type}.final_comp.vcf \
cat $(dirname {snakemake.output.done})/{snakemake.wildcards.me_type}.final_comp.vcf \
| bgzip -c \
> {snakemake.output.vcf}
tabix -f {snakemake.output.vcf}
"""
)
