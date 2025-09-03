# -*- coding: utf-8 -*-
"""Wrapper for running VCF2MAF incl VEP variant annotation"""

import os

from snakemake.shell import shell

args = getattr(snakemake.params, "args", {})

vcf_to_table = os.path.join(os.path.dirname(os.path.realpath(__file__)), "vcf_to_table.py")
if args["somatic_variant_annotation_tool"] == "vep":
    vcf_to_table_config = os.path.join(os.path.dirname(os.path.realpath(__file__)), "vep.yaml")
else:
    raise Exception(
        "vcf to maf conversion error: unimplemented conversion from annotation tool {}".format(
            args["somatic_variant_annotation_tool"]
        )
    )

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
    --NCBI_Build {args[ncbi_build]} --Center "{args[Center]}" \
    --vcf-tumor-id {args[tumor_sample]} --tumor-id {args[tumor_id]} \
    --vcf-normal-id {args[normal_sample]} --normal-id {args[normal_id]} \
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
