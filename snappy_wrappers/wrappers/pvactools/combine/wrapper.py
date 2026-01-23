# -*- coding: utf-8 -*-
"""
Wrapper to combine a somatic variant vcf with RNA expression values & pileup of RNA data at variant loci.
"""

import os

from snakemake.shell import shell

__author__ = "Eric Blanc"
__email__ = "eric.blanc@bih-charite.de"

args = getattr(snakemake.params, "args", {})

samples = args["tumor_sample"]
libraries = args["tumor_library"]
if "normal_sample" in args and "normal_library" in args:
    samples = args["normal_sample"] + " " + samples
    libraries = args["normal_library"] + "," + libraries

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
    --remove-duplicates --do-not-normalize --multiallelics -any \
    {snakemake.input.vcf} \
| bcftools view \
    --no-update \
    --samples "{libraries}" \
| bcftools reheader \
    --samples <(echo "{samples}" | tr ' ' '\n') \
| python {combine_script} \
    --sample {args[tumor_sample]} \
    {args[extra_args]} \
    --combine {snakemake.output.vcf}
tabix {snakemake.output.vcf}

pushd $(dirname {snakemake.output.vcf})
md5sum $(basename {snakemake.output.vcf}) >$(basename {snakemake.output.vcf}).md5
md5sum $(basename {snakemake.output.vcf_tbi}) >$(basename {snakemake.output.vcf_tbi}).md5
"""
)
