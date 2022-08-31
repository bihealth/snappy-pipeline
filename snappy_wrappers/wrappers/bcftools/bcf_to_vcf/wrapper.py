# -*- coding: utf-8 -*-
"""Wrapper for running bcftools convert - BCF to compressed VCF.
"""

from snakemake.shell import shell

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

bcftools convert \
    ---output-type z \
    --output {snakemake.output.vcf} \
    {snakemake.input.bcf}
tabix -f {snakemake.output.vcf}

pushd $(dirname {snakemake.output.vcf})
md5sum $(basename {snakemake.output.vcf}) > $(basename {snakemake.output.vcf_md5})
md5sum $(basename {snakemake.output.tbi}) > $(basename {snakemake.output.tbi_md5})
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
