from snappy_wrappers.snappy_wrapper import ShellWrapper

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

ShellWrapper(snakemake).run(
    r"""
set -x

sample_count=$(echo {snakemake.input.bcf} | tr ' ' '\n' | wc -l)

if [[ $sample_count -eq 1 ]]; then
    # If a single sample, there is no need to merge.
    bcftools view \
        -O z \
        -o {snakemake.output.vcf} \
        {snakemake.input.bcf}
else
    bcftools merge \
        -m id \
        -O z \
        -o {snakemake.output.vcf} \
        {snakemake.input.bcf}
fi
tabix -f {snakemake.output.vcf}
"""
)
