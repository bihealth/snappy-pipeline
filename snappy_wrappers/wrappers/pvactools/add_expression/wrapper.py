# -*- coding: utf-8 -*-
"""Wrapper to add RNA expression values & pileup of RNA data at variant loci to a somatic vcf"""

import os

from snappy_wrappers.snappy_wrapper import ShellWrapper

__author__ = "Eric Blanc"
__email__ = "eric.blanc@bih-charite.de"

args = getattr(snakemake.params, "args", {})
extra_args = " ".join(args.get("extra_args", []))

combine_script = os.path.join(os.path.dirname(__file__), "add_expression.py")

expression = []
if pileup := getattr(snakemake.input, "pileup", None):
    expression.append(f"--pileup {pileup}")
if gene_tpms := getattr(snakemake.input, "gene_tpms", None):
    expression.append(f"--gene-tpms {gene_tpms}")
if transcript_tpms := getattr(snakemake.input, "transcript_tpms", None):
    expression.append(f"--transcript-tpms {transcript_tpms}")
    if duplicates := getattr(snakemake.input, "duplicate_transcripts_table", None):
        expression.append(f"--duplicates {duplicates}")
if gene_tpms or transcript_tpms:
    expression.append(f"--format {args['format']}")
expression = " ".join(expression)

ShellWrapper(snakemake).run(
    r"""
bcftools norm \
    --rm-dup exact --do-not-normalize --multiallelics -any \
    {snakemake.input.annotated} \
| python {combine_script} \
    --sample {args[tumor_sample]} \
    {expression} \
    {extra_args} \
    --output {snakemake.output.vcf}
tabix {snakemake.output.vcf}
"""
)
