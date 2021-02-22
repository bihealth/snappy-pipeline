# -*- coding: utf-8 -*-
"""Wrapper for running Melt genotype
"""

from snakemake.shell import shell

__author__ = "Manuel Holtgrewe"
__email__ = "manuel.holtgrewe@bihealth.de"

shell(
    r"""
# -----------------------------------------------------------------------------
# Redirect stderr to log file by default and enable printing executed commands
exec 2> >(tee -a "{snakemake.log}")
set -x
# -----------------------------------------------------------------------------

melt_mei -Xmx2G Genotype \
    -h {snakemake.config[static_data_config][reference][path]} \
    -bamfile {snakemake.input.bam} \
    -p $(dirname {snakemake.input.done}) \
    -t $(melt_mei_path)/me_refs/{snakemake.config[step_config][wgs_mei_calling][melt][me_refs_infix]}/{snakemake.wildcards.me_type}_MELT.zip \
    -w $(dirname {snakemake.output.done})

"""
)
