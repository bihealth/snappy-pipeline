"""WhatsHap read-based phasing of the one sample of a VCF with its alignments

The phase goes into the ``HP`` tag of GATK 3 ReadBackedPhasing (``--tag=HP``), which pVACseq reads
for the proximal variants of neoepitope prediction; the genotypes stay unchanged."""

from snappy_wrappers.snappy_wrapper import ShellWrapper

ShellWrapper(snakemake).run(
    r"""
sample=$(bcftools query --list-samples {snakemake.input.vcf} | head -n 1)

whatshap phase \
    --tag=HP \
    {extra_arguments} \
    --reference {snakemake.input.reference} \
    --ignore-read-groups --sample "$sample" \
    --output $TMPDIR/phased.vcf \
    {snakemake.input.vcf} \
    {snakemake.input.bam}

bgzip -c $TMPDIR/phased.vcf > {snakemake.output.vcf}
tabix -f {snakemake.output.vcf}
"""
)
