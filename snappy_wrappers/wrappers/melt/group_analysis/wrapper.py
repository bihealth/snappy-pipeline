# -*- coding: utf-8 -*-
"""Wrapper for running Melt group_analysis
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

# TODO: read length should not be hard-coded but data set property
melt_mei -Xmx2G GroupAnalysis \
    -h {snakemake.config[static_data_config][reference][path]} \
    -t $(melt_mei_path)/me_refs/{snakemake.config[step_config][wgs_mei_calling][melt][me_refs_infix]}/{snakemake.wildcards.me_type}_MELT.zip \
    -n $(melt_mei_path)/{snakemake.config[step_config][wgs_mei_calling][melt][genes_file]} \
    $(if echo "{snakemake.config[step_config][wgs_mei_calling][melt][me_refs_infix]}" | grep -i hg19; then
        echo -v $(melt_mei_path)/prior_files/{snakemake.wildcards.me_type}.1KGP.sites.vcf;
    fi) \
    -w $(dirname {snakemake.output.done}) \
    -r 150 \
    -discoverydir $(dirname $(echo {snakemake.input} | cut -d ' ' -f 1))
"""
)
