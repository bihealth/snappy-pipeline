# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for bcftools roh
"""

from snakemake.shell import shell

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>"

shell(
    r"""
set -x

export TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT

# Also pipe stderr to log file
if [[ -n "{snakemake.log}" ]]; then
    if [[ "$(set +e; tty; set -e)" != "" ]]; then
        rm -f "{snakemake.log}" && mkdir -p $(dirname {snakemake.log})
        exec 2> >(tee -a "{snakemake.log}" >&2)
    else
        rm -f "{snakemake.log}" && mkdir -p $(dirname {snakemake.log})
        echo "No tty, logging disabled" >"{snakemake.log}"
    fi
fi

# Perform ROH Calling -----------------------------------------------------------------------------

bcftools roh \
    $(if [[ "{snakemake.config[step_config][roh_calling][bcftools_roh][af_tag]}" != "None" ]]; then
        echo --AF-tag "{snakemake.config[step_config][roh_calling][bcftools_roh][af_tag]}"
    fi) \
    $(if [[ "{snakemake.config[step_config][roh_calling][bcftools_roh][gts_only]}" != "None" ]]; then
        echo --GTs-only "{snakemake.config[step_config][roh_calling][bcftools_roh][gts_only]}"
    fi) \
    $(if [[ "{snakemake.config[step_config][roh_calling][bcftools_roh][ignore_homref]}" != "False" ]]; then
        echo --ignore-homref
    fi) \
    $(if [[ "{snakemake.config[step_config][roh_calling][bcftools_roh][skip_indels]}" != "False" ]]; then
        echo --skip-indels
    fi) \
    $(if [[ "{snakemake.config[step_config][roh_calling][bcftools_roh][rec_rate]}" != "None" ]]; then
        echo --rec-rate "{snakemake.config[step_config][roh_calling][bcftools_roh][rec_rate]}"
    fi) \
    -o {snakemake.output.txt} \
    -O srz \
    --threads $(({snakemake.config[step_config][roh_calling][bcftools_roh][threads]} - 1)) \
    {snakemake.input.vcf}

# Generate MD5 sum files

pushd $(dirname {snakemake.output.txt})
md5sum $(basename {snakemake.output.txt}) >$(basename {snakemake.output.txt}).md5
"""
)
