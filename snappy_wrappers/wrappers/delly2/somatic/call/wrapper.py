# -*- coding: utf-8 -*-
"""Wrapper for running Delly2's call step on tumor/matched-normal pairs"""

from snakemake.shell import shell

__author__ = "Nina Thiessen"
__email__ = "nina.thiessen@bih-charite.de"

exclude_str = ""
s = snakemake.config["step_config"]["somatic_wgs_sv_calling"]["delly2"]["path_exclude_tsv"]
if s is not None:
    exclude_str = "--exclude {}".format(s)

shell(
    r"""
# -----------------------------------------------------------------------------
# Redirect stderr to log file by default and enable printing executed commands
exec &> >(tee -a "{snakemake.log}")
set -x
# -----------------------------------------------------------------------------

delly call \
    {exclude_str} \
    --map-qual 1 \
    --qual-tra 20 \
    --genome {snakemake.config[static_data_config][reference][path]} \
    --outfile {snakemake.output.bcf} \
    {snakemake.input.tumor_bam} \
    {snakemake.input.normal_bam}

tabix -f {snakemake.output.bcf}

pushd $(dirname {snakemake.output.bcf})
md5sum $(basename {snakemake.output.bcf}) >$(basename {snakemake.output.bcf}).md5
md5sum $(basename {snakemake.output.bcf}).csi >$(basename {snakemake.output.bcf}).csi.md5
"""
)
