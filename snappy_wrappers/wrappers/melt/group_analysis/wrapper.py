from math import ceil

from snappy_wrappers.snappy_wrapper import ShellWrapper

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

args = getattr(snakemake.params, "args", {})

reference = snakemake.input.reference

# As an estimate, give 75% of total memory as usable to JVM
mem_jvm = ceil(snakemake.resources.mem_mb * 0.75)

ShellWrapper(snakemake).run(
    r"""

JAR={snakemake.input.jar}
ME_REFS={snakemake.input.me_refs}
ME_INFIX={args[me_refs_infix]}

# GroupAnalysis reads the IndivAnalysis output of all pedigree members from one directory,
# but each library's IndivAnalysis wrote to its own: link them together (the .done files are hidden)
discovery=$TMPDIR/discovery
mkdir -p $discovery
for done in {snakemake.input.indiv_analysis}; do
    ln -s $(realpath $(dirname $done))/* $discovery/
done

java -Xmx{mem_jvm}m -jar $JAR \
    GroupAnalysis \
    -h {reference} \
    -t $ME_REFS/$ME_INFIX/{args[me_type]}_MELT.zip \
    $(if [[ $ME_REFS == *37* ]] || [[ $ME_REFS == *hg19* ]]; then
        echo -v $ME_REFS/../../prior_files/{args[me_type]}.1KGP.sites.vcf;
    fi) \
    -w $(dirname {snakemake.output.done}) \
    -r {args[read_length]} \
    -n {snakemake.input.genes} \
    -discoverydir $discovery
"""
)
