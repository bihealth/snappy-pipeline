# -*- coding: utf-8 -*-
"""Wrapper for running VCF2MAF incl VEP variant annotation"""

import os
from typing import TYPE_CHECKING

from snappy_wrappers.snappy_wrapper import ShellWrapper

if TYPE_CHECKING:
    from snakemake.iocontainers import snakemake

args = getattr(snakemake.params, "args", {})

samples = f"--vcf-tumor-id {args['tumor_sample']} --tumor-id {args['tumor_id']}"
if args["normal_sample"] and args["normal_id"]:
    samples += f" --vcf-normal-id {args['normal_sample']} --normal-id {args['normal_id']}"

vcf_to_table = os.path.join(os.path.dirname(os.path.realpath(__file__)), "vcf_to_table.py")
config_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), args["somatic_variant_annotation_tool"] + ".yaml")
if os.path.exists(config_path) and os.path.isfile(config_path):
    vcf_to_table_config = config_path
else:
    raise Exception(
        "vcf to maf conversion error: unimplemented conversion from annotation tool {}".format(
            args["somatic_variant_annotation_tool"]
        )
    )

ShellWrapper(snakemake).run(
    r"""
python {vcf_to_table} \
    --config {vcf_to_table_config} \
    --debug --unique --title \
    --NCBI_Build {args[ncbi_build]} --Center "{args[Center]}" \
    {samples} \
    {snakemake.input.vcf} {snakemake.output.maf}
"""
)
