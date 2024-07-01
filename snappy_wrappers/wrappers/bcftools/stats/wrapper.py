# -*- coding: utf-8 -*-
"""Wrapper for running Canvas in somatic variant calling mode on WGS data"""

from snakemake.shell import shell

__author__ = "Manuel Holtgrewe"
__email__ = "manuel.holtgrewe@bih-charite.de"

shell(
    r"""
# -----------------------------------------------------------------------------
# Redirect stderr to log file by default and enable printing executed commands
exec 2> >(tee -a "{snakemake.log}")
set -x
# -----------------------------------------------------------------------------

module purge
module load HTSlib/1.3.1-foss-2015a
module load BCFtools/1.3.1-foss-2015a

bcftools stats {snakemake.input} \
> {snakemake.output.txt}

pushd $(dirname {snakemake.output.txt})
md5sum $(basename {snakemake.output.txt}) >$(basename {snakemake.output.txt}).md5
"""
)
