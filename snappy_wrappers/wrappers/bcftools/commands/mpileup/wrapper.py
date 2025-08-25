# -*- coding: utf-8 -*-
"""
Wrapper for finding variants using bcftools mpileup

It is meant to be used in conjunction with other bcftools commands, such as call & filter

Mandatory snakemake.input: bams
Optional snakemake.input: fasta_ref, read_groups, regions_file, samples_file, targets_file

Mandatory snakemake.params.args: extra_args
Optional snakemake.params.args: index

Mandatory snakemake.output: vcf
Optional snakemake.output:
"""

__author__ = "Eric Blanc <eric.blanc@bih-charite.de>"

import os
import sys

# The following is required for being able to import snappy_wrappers modules
# inside wrappers. When the wrappers have their own python environment, messing
# with the path is necessary.
base_dir = os.path.normpath(os.path.dirname(__file__))
while os.path.basename(base_dir) != "snappy_wrappers":
    base_dir = os.path.dirname(base_dir)
sys.path.insert(0, os.path.dirname(base_dir))

from snappy_wrappers.snappy_wrapper import ShellWrapper, BcftoolsWrapper, BcftoolsCommand

args = getattr(snakemake.params, "args", {})

extra_files = []
for additional_input in BcftoolsWrapper.additional_input_files[BcftoolsCommand.MPILEUP]:
    filename = getattr(snakemake.input, additional_input, None)
    if filename is not None:
        extra_files.append(
            f"--{additional_input.replace('_', '-')} {filename}"
        )

cmd = r"""
bcftools mpileup \
    {extra_args} \
    {extra_files} \
    --output {snakemake.output.vcf} \
    {snakemake.input.bams}
""".format(
    snakemake=snakemake,
    extra_args=" \\\n    ".join(args["extra_args"]) if args.get("extra_args", None) is not None else "",
    extra_files=" \\\n    ".join(extra_files) if extra_files else "",
)

if args.get("index", False):
    cmd += "\n" + BcftoolsWrapper.index_command.format(snakemake=snakemake)

ShellWrapper(snakemake).run(cmd)
