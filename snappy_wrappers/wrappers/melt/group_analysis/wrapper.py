from snappy_wrappers.snappy_wrapper import ShellWrapper

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

args = getattr(snakemake.params, "args", {})

reference = snakemake.input.reference
individual = snakemake.input.indiv_analysis[0]

ShellWrapper(snakemake).run(
    r"""

JAR={snakemake.input.jar}
ME_REFS={snakemake.input.me_refs}
ME_INFIX={args[me_refs_infix]}

java -jar -Xmx13G -jar $JAR \
    GroupAnalysis \
    -h {reference} \
    -t $ME_REFS/$ME_INFIX/{args[me_type]}_MELT.zip \
    $(if [[ $ME_REFS == *37* ]] || [[ $ME_REFS == *hg19* ]]; then
        echo -v $ME_REFS/../../prior_files/{args[me_type]}.1KGP.sites.vcf;
    fi) \
    -w $(dirname {snakemake.output.done}) \
    -r 150 \
    -n {snakemake.input.genes} \
    -discoverydir $(dirname {individual})
"""
)
