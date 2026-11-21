from snappy_wrappers.snappy_wrapper import ShellWrapper
from snappy_wrappers.tools.java import java_heap_mb

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

args = getattr(snakemake.params, "args", {})

reference = snakemake.input.reference

mem_jvm = java_heap_mb(snakemake)

ShellWrapper(snakemake).run(
    r"""
JAR={snakemake.input.jar}
ME_REFS={snakemake.input.me_refs}
ME_INFIX={args[me_refs_infix]}

# GroupAnalysis reads the IndivAnalysis output of all pedigree members from one directory,
# but each library's IndivAnalysis wrote to its own: link their files together (the glob leaves
# out Snakemake's hidden .snakemake_timestamp files)
discovery=$TMPDIR/discovery
mkdir -p $discovery
for dir in {snakemake.input.indiv_analysis}; do
    ln -s $(realpath $dir)/* $discovery/
done

java -Xmx{mem_jvm}m -jar $JAR \
    GroupAnalysis \
    -h {reference} \
    -t $ME_REFS/$ME_INFIX/{args[me_type]}_MELT.zip \
    $(if [[ $ME_REFS == *37* ]] || [[ $ME_REFS == *hg19* ]]; then
        echo -v $ME_REFS/../../prior_files/{args[me_type]}.1KGP.sites.vcf
    fi) \
    -w $(dirname {snakemake.output.pre_geno}) \
    -r {args[read_length]} \
    -n {snakemake.input.genes} \
    -discoverydir $discovery
"""
)
