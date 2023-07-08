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
    --exom-bedfile {snakemake.config[step_config][somatic_variant_checking][target_regions]} \
    --padding {snakemake.config[step_config][somatic_variant_checking][padding]} \
    --minimal {snakemake.config[step_config][somatic_variant_checking][minimal_support_read]} \
    --limited {snakemake.config[step_config][somatic_variant_checking][limited_support_read]} \
    --ignore-regions {snakemake.config[step_config][somatic_variant_checking][ignore_regions]} \
    --AF {snakemake.config[step_config][somatic_variant_checking][AF_ID]} \
    --output {snakemake.output.json}

pushd $(dirname {snakemake.output.json})
md5sum $(basename {snakemake.output.json}) > $(basename {snakemake.output.json_md5})
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
