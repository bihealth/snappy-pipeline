# -*- coding: utf-8 -*-
"""Wrapper for pedigree extraction from popdel cohort-wide call step."""

from snakemake.shell import shell

__author__ = "Manuel Holtgrewe"
__email__ = "manuel.holtgrewe@bih-charite.de"

shell(
    r"""
# -----------------------------------------------------------------------------
# Redirect stderr to log file by default and enable printing executed commands
exec &> >(tee -a "{snakemake.log}")
set -x
# -----------------------------------------------------------------------------

export TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT

echo '{snakemake.params.ped_members}' \
| tr ' ' '\n' \
>$TMPDIR/samples.txt

bcftools view \
    --samples-file $TMPDIR/samples.txt \
    --output-type u \
    {snakemake.input.vcf} \
| bcftools view \
    --output-file {snakemake.output.vcf} \
    --output-type z \
    --include '(GT !~ "\.") && (GT ~ "1")'

tabix -f {snakemake.output.vcf}

pushd $(dirname {snakemake.output.vcf})
md5sum $(basename {snakemake.output.vcf}) >$(basename {snakemake.output.vcf}).md5
md5sum $(basename {snakemake.output.vcf}).tbi >$(basename {snakemake.output.vcf}).tbi.md5
"""
)
