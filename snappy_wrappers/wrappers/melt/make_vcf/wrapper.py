from math import ceil

from snappy_wrappers.snappy_wrapper import ShellWrapper

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

args = getattr(snakemake.params, "args", {})

# As an estimate, give 75% of total memory as usable to JVM
mem_jvm = ceil(snakemake.resources.mem_mb * 0.75)

ShellWrapper(snakemake).run(
    r"""
JAR={snakemake.input.jar}
ME_REFS={snakemake.input.me_refs}
ME_INFIX={args[me_refs_infix]}

out=$(dirname {snakemake.output.vcf})
genotype_dir=$(dirname {snakemake.input.genotype} | head -n 1)
ls $genotype_dir/*.{args[me_type]}.tsv \
    | sort \
        > {snakemake.output.list_txt}

# TODO: allowing 100% no-call wise?
java -Xmx{mem_jvm}m -jar $JAR \
    MakeVCF \
    -genotypingdir $genotype_dir \
    -h {snakemake.input.reference} \
    -j {args[max_no_call_percent]} \
    -t $ME_REFS/$ME_INFIX/{args[me_type]}_MELT.zip \
    -p $(dirname {snakemake.input.pre_geno}) \
    -w $out \
    -o $out

## Make file more VCF conforming XXX TODO XXX
#perl -p -e 's/ID=GL,Number=3/ID=GL,Number=G/' $out/{args[me_type]}.final_comp.vcf \
cat $out/{args[me_type]}.final_comp.vcf \
    | bgzip -c \
        > {snakemake.output.vcf}
tabix -f {snakemake.output.vcf}
"""
)
