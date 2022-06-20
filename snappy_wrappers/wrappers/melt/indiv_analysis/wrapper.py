# -*- coding: utf-8 -*-
"""Wrapper for running Melt IndivAnalysis
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

export MALLOC_ARENA_MAX=4

melt_mei -Xmx18G IndivAnalysis \
    -b hs37d5/NC_007605 \
    -c 30 \
    -h {snakemake.config[static_data_config][reference][path]} \
    -t $(melt_mei_path)/me_refs/{snakemake.config[step_config][wgs_mei_calling][melt][me_refs_infix]}/{snakemake.wildcards.me_type}_MELT.zip \
    -w $(dirname {snakemake.output.done}) \
    -r 150 \
    -bamfile {snakemake.input.orig_bam}
"""
)
