# -*- coding: utf-8 -*-
"""Wrapper for summarize information of vcf file after
somatic variant calling step
"""
import os

from snakemake.shell import shell

__author__ = "Pham Gia Cuong"
__email__ = "pham.gia-cuong@bih-charite.de"

base_dir = os.path.dirname(os.path.realpath(__file__))
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

python {base_dir}/summarize-vcf.py --rawvcf {snakemake.input.full_vcf} \
    --filtered-vcf {snakemake.input.passed_vcf} \
    --sample {snakemake.wildcards.tumor_library} \
    --exom-bedfile {snakemake.config[step_config][somatic_variant_checking][target_regions]} \
    --padding {snakemake.config[step_config][somatic_variant_checking][padding]} \
    --minimal {snakemake.config[step_config][somatic_variant_checking][minimal_support_read]} \
    --limited {snakemake.config[step_config][somatic_variant_checking][limited_support_read]} \
    --ignore-regions {snakemake.config[step_config][somatic_variant_checking][ignore_regions]} \
    --variant-allele-frequency-id {snakemake.config[step_config][somatic_variant_checking][variant_allele_frequency_id]} \
    --output {snakemake.output.json}

pushd $(dirname {snakemake.output.json})
md5sum $(basename {snakemake.output.json}) > $(basename {snakemake.output.json_md5})
popd
    """
)

# Compute MD5 sums of logs
shell(
    r"""
pushd $(dirname {snakemake.log.log})
md5sum $(basename {snakemake.log.log}) > $(basename {snakemake.log.log_md5})
md5sum $(basename {snakemake.log.conda_list}) > $(basename {snakemake.log.conda_list_md5})
md5sum $(basename {snakemake.log.conda_info}) > $(basename {snakemake.log.conda_info_md5})
"""
)
