# -*- coding: utf-8 -*-
"""Wrapper for running VCF2MAF incl VEP variant annotation
"""

import os

from snakemake.shell import shell

step = snakemake.config["pipeline_step"]["name"]
config = snakemake.config["step_config"][step]

vcf_to_table = os.path.join(os.path.dirname(os.path.realpath(__file__)), "vcf_to_table.py")
vcf_to_table_config = os.path.join(os.path.dirname(os.path.realpath(__file__)), "config.yaml")

params = snakemake.params.args

shell(
    r"""
set -x

# Also pipe everything to log file
if [[ -n "{snakemake.log.log}" ]]; then
    if [[ "$(set +e; tty; set -e)" != "" ]]; then
        rm -f "{snakemake.log.log}" && mkdir -p $(dirname {snakemake.log.log})
        exec &> >(tee -a "{snakemake.log.log}" >&2)
    else
        rm -f "{snakemake.log.log}" && mkdir -p $(dirname {snakemake.log.log})
        echo "No tty, logging disabled" >"{snakemake.log.log}"
    fi
fi

# Write out information about conda installation.
conda list >{snakemake.log.conda_list}
conda info >{snakemake.log.conda_info}
md5sum {snakemake.log.conda_list} >{snakemake.log.conda_list_md5}
md5sum {snakemake.log.conda_info} >{snakemake.log.conda_info_md5}

python {vcf_to_table} \
    --config {vcf_to_table_config} \
    --debug --unique --title \
    --NCBI_Build {config[vcf2maf][ncbi_build]} --Center {config[vcf2maf][Center]} \
    --vcf-tumor-id {params[tumor_sample]} --tumor-id {params[tumor_id]} \
    --vcf-normal-id {params[normal_sample]} --normal-id {params[normal_id]} \
    {snakemake.input.vcf} {snakemake.output.maf}

pushd $(dirname {snakemake.output.maf})
md5sum $(basename {snakemake.output.maf}) > $(basename {snakemake.output.maf}).md5
popd
"""
)

# Compute MD5 sums of logs.
shell(
    r"""
md5sum {snakemake.log.log} >{snakemake.log.log_md5}
"""
)
