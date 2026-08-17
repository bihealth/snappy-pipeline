from snappy_wrappers.snappy_wrapper import ShellWrapper

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

ShellWrapper(snakemake).run(
    r"""
# Filter input VCF files to those with any lines
input=
for vcf in {snakemake.input.vcf}; do
    lines=$(zcat $vcf | wc -l)
    if [[ $lines -ne 0 ]]; then
        input="$input $vcf"
    fi
done

bcftools concat -a $input \
| bgzip -c \
> {snakemake.output.vcf}
tabix -f {snakemake.output.vcf}
"""
)
