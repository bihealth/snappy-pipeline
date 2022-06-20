# -*- coding: utf-8 -*-
"""Wrapper for running ERDS+SV2 merge genotypes step
"""

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

bcftools merge \
    -m id \
    -O z \
    -o {snakemake.output.vcf} \
    {snakemake.input}

$(which tabix) --version

$(which tabix) -f {snakemake.output.vcf}

pushd $(dirname {snakemake.output.vcf})
md5sum $(basename {snakemake.output.vcf}) >$(basename {snakemake.output.vcf}).md5
md5sum $(basename {snakemake.output.vcf}).tbi >$(basename {snakemake.output.vcf}).tbi.md5
"""
)
