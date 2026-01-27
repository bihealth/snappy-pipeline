# -*- coding: utf-8 -*-
"""
Wrapper to combine a somatic variant vcf with RNA expression values & pileup of RNA data at variant loci.
"""

import os

from snakemake.shell import shell

__author__ = "Eric Blanc"
__email__ = "eric.blanc@bih-charite.de"

args = getattr(snakemake.params, "args", {})

combine_script = os.path.join(os.path.dirname(__file__), "combine.py")

# Renaming libraries to samples, and ensure normal before tumor
# (bcftools view ensures order of libraries, bcftools reheader renames)
samples = args["tumor_sample"]
libraries = args["tumor_library"]
if "normal_sample" in args and "normal_library" in args:
    samples = args["normal_sample"] + " " + samples
    libraries = args["normal_library"] + "," + libraries

expression = []
if pileup := getattr(snakemake.input, "pileup", None):
    expression.append(f"--pileup {pileup}")
if gene_tpms := getattr(snakemake.input, "gene_tpms", None):
    expression.append(f"--gene-tpms {gene_tpms}")
if transcript_tpms := getattr(snakemake.input, "transcript_tpms", None):
    expression.append(f"--transcript-tpms {transcript_tpms}")
expression = " ".join(expression)

shell(
    r"""
# -----------------------------------------------------------------------------
# Redirect stderr to log file by default and enable printing executed commands
# Also pipe stderr to log file
if [[ -n "{snakemake.log.log}" ]]; then
    if [[ "$(set +e; tty; set -e)" != "" ]]; then
        rm -f "{snakemake.log.log}" && mkdir -p $(dirname {snakemake.log.log})
        exec 2> >(tee -a "{snakemake.log.log}" >&2)
    else
        rm -f "{snakemake.log.log}" && mkdir -p $(dirname {snakemake.log.log})
        echo "No tty, logging disabled" >"{snakemake.log.log}"
    fi
fi
# -----------------------------------------------------------------------------
# Create auto-cleaned temporary directory
export TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT

set -x
# -----------------------------------------------------------------------------

# Write out information about conda installation
conda list > {snakemake.log.conda_list}
conda info > {snakemake.log.conda_info}

bcftools norm \
    --rm-dup exact --do-not-normalize --multiallelics -any \
    {snakemake.input.annotated} \
| bcftools view \
    --no-update \
    --samples "{libraries}" \
| bcftools reheader \
    --samples <(echo "{samples}" | tr ' ' '\n') \
| python {combine_script} \
    --sample {args[tumor_sample]} \
    {expression} \
    {args[extra_args]} \
    --output {snakemake.output.vcf}
tabix {snakemake.output.vcf}

pushd $(dirname {snakemake.output.vcf})
md5sum $(basename {snakemake.output.vcf}) >$(basename {snakemake.output.vcf}).md5
md5sum $(basename {snakemake.output.vcf}).tbi >$(basename {snakemake.output.vcf}).tbi.md5
"""
)
