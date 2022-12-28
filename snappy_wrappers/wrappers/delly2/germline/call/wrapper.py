# -*- coding: utf-8 -*-
"""Wrapper for running Delly2's call step
"""

from snakemake.shell import shell

__author__ = "Manuel Holtgrewe"
__email__ = "manuel.holtgrewe@bih-charite.de"

exclude_str = ""
if "wgs_sv_calling" in snakemake.config["step_config"]:
    s = snakemake.config["step_config"]["wgs_sv_calling"]["delly2"]["path_exclude_tsv"]
else:
    s = snakemake.config["step_config"]["sv_calling_targeted"]["delly2"]["path_exclude_tsv"]
if s:
    exclude_str = f"--exclude {s}"

shell(
    r"""
# -----------------------------------------------------------------------------
# Redirect stderr to log file by default and enable printing executed commands
exec &> >(tee -a "{snakemake.log}")
set -x
# -----------------------------------------------------------------------------

delly call \
    --map-qual {snakemake.config[step_config][wgs_sv_calling][delly2][map_qual]} \
    --qual-tra {snakemake.config[step_config][wgs_sv_calling][delly2][qual_tra]} \
    --geno-qual {snakemake.config[step_config][wgs_sv_calling][delly2][geno_qual]} \
    --mad-cutoff {snakemake.config[step_config][wgs_sv_calling][delly2][mad_cutoff]} \
    --genome {snakemake.config[static_data_config][reference][path]} \
    --outfile {snakemake.output.bcf} \
    {exclude_str} \
    {snakemake.input.bam}

tabix -f {snakemake.output.bcf}

pushd $(dirname {snakemake.output.bcf})
md5sum $(basename {snakemake.output.bcf}) >$(basename {snakemake.output.bcf}).md5
md5sum $(basename {snakemake.output.bcf}).csi >$(basename {snakemake.output.bcf}).csi.md5
"""
)
