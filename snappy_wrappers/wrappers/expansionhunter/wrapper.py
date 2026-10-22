# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for ExpansionHunter: Snakemake wrapper.py"""

from typing import TYPE_CHECKING

from snappy_wrappers.snappy_wrapper import ShellWrapper

if TYPE_CHECKING:
    from snakemake.iocontainers import snakemake

args = getattr(snakemake.params, "args", {})

# Define argument sex if any (otherwise: female [default])
sex_argument = ""
valid_sex_list = ["female", "male"]
if args["sex"] in valid_sex_list:
    sex_argument = "--sex " + args["sex"]


ShellWrapper(snakemake).run(
    r"""
# ExpansionHunter also writes <prefix>_realigned.bam, which is not declared
prefix=$TMPDIR/expansionhunter

ExpansionHunter --reads {snakemake.input.bam} \
        --reference {snakemake.input.reference} \
        --variant-catalog {snakemake.input.repeat_catalog} \
        --output-prefix $prefix {sex_argument} {extra_arguments}

cp $prefix.json {snakemake.output.json}
cp $prefix.vcf {snakemake.output.vcf}
"""
)
