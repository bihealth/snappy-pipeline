# -*- coding: utf-8 -*-
"""Wrapper for running Melt preprocessing
"""

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

ln -sr {snakemake.input.bam} {snakemake.output.orig_bam}
ln -sr {snakemake.input.bai} {snakemake.output.orig_bai}

melt_mei -Xmx2G Preprocess \
    Preprocess \
    -bamfile {snakemake.output.orig_bam} \
    -h {snakemake.config[static_data_config][reference][path]}
"""
)
