# -*- coding: utf-8 -*-
"""Wrapper code for GATK PhaseByTransmission"""

from math import ceil
from typing import TYPE_CHECKING

from snappy_wrappers.snappy_wrapper import ShellWrapper

if TYPE_CHECKING:
    from snakemake.iocontainers import snakemake

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

args = getattr(snakemake.params, "args", {})

java_mem_mb = ceil(snakemake.resources.mem_mb * 0.75)
father_allele_first = "--FatherAlleleFirst" if args["father_allele_first"] else ""

ShellWrapper(snakemake).run(
    r"""
which tabix
which bcftools
which gatk_nonfree

# Extract pedigree only
members=$(cut -f 2 {snakemake.input.ped} | tr '\n' ',' | sed -e 's/,$//'g)
bcftools view \
    --threads {snakemake.threads} \
    -s "$members" \
    {snakemake.input.vcf} \
    -O z \
    -o $TMPDIR/trio_only.vcf.gz
tabix -f $TMPDIR/trio_only.vcf.gz

# Call wrapper
MALLOC_ARENA_MAX=4 \
    gatk_nonfree \
    -Xmx{java_mem_mb}m \
    -Djava.io.tmpdir=$TMPDIR \
    --analysis_type PhaseByTransmission \
    -nct 1 \
    --pedigreeValidationType {args[pedigree_validation_type]} \
    {father_allele_first} \
    --DeNovoPrior {args[de_novo_prior]} \
    {extra_arguments} \
    --pedigree {snakemake.input.ped} \
    --variant $TMPDIR/trio_only.vcf.gz \
    --out {snakemake.output.vcf} \
    --reference_sequence {snakemake.input.reference}

tabix -f {snakemake.output.vcf}
"""
)
