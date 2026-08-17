from snappy_wrappers.snappy_wrapper import ShellWrapper

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

args = getattr(snakemake.params, "args", {})

ShellWrapper(snakemake).run(
    r"""
bcftools stats \
    --split-by-ID \
    -s {args[donor_library_name]} \
    {snakemake.input.vcf} \
> {snakemake.output.txt}
"""
)
