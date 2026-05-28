# -*- coding: utf-8 -*-
"""Wrapper for running bcftools convert - gVCF to VCF."""

from typing import TYPE_CHECKING

from snappy_wrappers.snappy_wrapper import ShellWrapper

if TYPE_CHECKING:
    from snakemake.iocontainers import snakemake

args = getattr(snakemake.params, "args", {})
reference_path = args["reference_path"]

ShellWrapper(snakemake).run(
    r"""
# Method checks if VCF contains sample
check_vcf() {{
    # Variables
    # vcf=$1
    # sample=$2

    # Check
    if bcftools query --list-samples $1 | grep --quiet --word-regexp $2; then
       return 0
    else
        echo "VCF header doesn't contain sample '$2': $1"
        echo "Samples:" $(bcftools query --list-samples $1)
        exit 1
    fi
}}

# Convert gVCF to VCF, filter at least one allele
bcftools convert --gvcf2vcf \
        --output-type u \
        --fasta-ref {reference_path} \
        {snakemake.params.args[input]} \
| bcftools view --no-update --min-ac 1 \
        --output-type z \
        --output {snakemake.output.vcf}
tabix -f {snakemake.output.vcf}

# Validate VCF: contains all expected samples
while read sample; do
    check_vcf {snakemake.output.vcf} $sample
done < <(echo {snakemake.params.args[sample_names]})

"""
)
