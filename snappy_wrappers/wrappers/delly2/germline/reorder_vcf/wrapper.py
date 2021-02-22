# -*- coding: utf-8 -*-
"""Wrapper for running Delly2's calls tep
"""

from snakemake.shell import shell

__author__ = "Manuel Holtgrewe"
__email__ = "manuel.holtgrewe@bihealth.de"

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
    {snakemake.input.bcf} \
| bcftools view \
    --output-file {snakemake.output.vcf} \
    --output-type z \
    --include '(GT !~ "\.") && (GT ~ "1")'

tabix -s1 -b2 -e2 -f {snakemake.output.vcf}

pushd $(dirname {snakemake.output.vcf})
md5sum $(basename {snakemake.output.vcf}) >$(basename {snakemake.output.vcf}).md5
md5sum $(basename {snakemake.output.vcf}).tbi >$(basename {snakemake.output.vcf}).tbi.md5
"""
)
