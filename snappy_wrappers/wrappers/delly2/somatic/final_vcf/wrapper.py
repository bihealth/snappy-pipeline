# -*- coding: utf-8 -*-
"""Wrapper for running Delly2's "final vcf" step
"""

from snakemake.shell import shell

__author__ = "Nina Thiessen"
__email__ = "nina.thiessen@bih-charite.de"

shell(
    r"""
# -----------------------------------------------------------------------------
# Redirect stderr to log file by default and enable printing executed commands
exec &> >(tee -a "{snakemake.log}")
set -x
# -----------------------------------------------------------------------------

bcftools view
    -O z
    --outfile {snakemake.output.vcf} \
    {snakemake.input.bcf}

tabix -f {snakemake.output.vcf}

pushd $(dirname {snakemake.output.vcf})
md5sum $(basename {snakemake.output.vcf}) >$(basename {snakemake.output.vcf}).md5
md5sum $(basename {snakemake.output.vcf}).tbi >$(basename {snakemake.output.vcf}).tbi.md5
"""
)
