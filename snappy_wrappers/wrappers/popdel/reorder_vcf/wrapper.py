from snappy_wrappers.snappy_wrapper import ShellWrapper

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

args = getattr(snakemake.params, "args", {})

ShellWrapper(snakemake).run(
    r"""
echo '{args[ped_members]}' \
| tr ' ' '\n' \
>$TMPDIR/samples.txt

bcftools view \
    --samples-file $TMPDIR/samples.txt \
    --output-type u \
    {snakemake.input.vcf} \
| bcftools view \
    --output-file {snakemake.output.vcf} \
    --output-type z \
    --include '(GT !~ "\.") && (GT ~ "1")'

tabix -f {snakemake.output.vcf}
"""
)
