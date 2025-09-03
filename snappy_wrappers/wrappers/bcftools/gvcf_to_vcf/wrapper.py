# -*- coding: utf-8 -*-
"""Wrapper for running bcftools convert - gVCF to VCF."""

from typing import TYPE_CHECKING

from snakemake.shell import shell

if TYPE_CHECKING:
    from snakemake.script import snakemake

args = getattr(snakemake.params, "args", {})
reference_path = args["reference_path"]

shell(
    r"""
# -----------------------------------------------------------------------------
# Redirect stderr to log file by default and enable printing executed commands
exec 2> >(tee -a "{snakemake.log.log}")
set -x
# -----------------------------------------------------------------------------

# Write out information about conda installation
conda list > {snakemake.log.conda_list}
conda info > {snakemake.log.conda_info}

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

pushd $(dirname {snakemake.output.vcf})
md5sum $(basename {snakemake.output.vcf}) > $(basename {snakemake.output.vcf_md5})
md5sum $(basename {snakemake.output.vcf_tbi}) > $(basename {snakemake.output.vcf_tbi_md5})
"""
)

# Compute MD5 sums of logs
shell(
    r"""
md5sum {snakemake.log.log} > {snakemake.log.log_md5}
md5sum {snakemake.log.conda_list} > {snakemake.log.conda_list_md5}
md5sum {snakemake.log.conda_info} > {snakemake.log.conda_info_md5}
"""
)
