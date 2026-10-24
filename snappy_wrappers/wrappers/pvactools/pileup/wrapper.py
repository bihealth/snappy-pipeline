# -*- coding: utf-8 -*-
"""Wrapper to pileup RNA data at variant loci"""

from snappy_wrappers.snappy_wrapper import ShellWrapper

__author__ = "Eric Blanc"
__email__ = "eric.blanc@bih-charite.de"

args = getattr(snakemake.params, "args", {})

ShellWrapper(snakemake).run(
    r"""
export TMPDIR=$(realpath $TMPDIR)
trap "rm -rf $TMPDIR" EXIT

bcftools query \
    --format '%CHROM\t%POS0\t%END' \
    {snakemake.input.loci} \
    | bedtools merge -i - \
        > $TMPDIR/regions.bed

bcftools mpileup \
    --fasta-ref {snakemake.input.reference} \
    --regions-file $TMPDIR/regions.bed \
    {extra_arguments} \
    {snakemake.input.bam} \
    | bcftools reheader --samples <(echo "{args[tumor_sample]}") \
    | bcftools norm --multiallelics -both --remove-duplicates \
    | bcftools view -AA --output-type z --output {snakemake.output.vcf} --write-index=tbi
"""
)
