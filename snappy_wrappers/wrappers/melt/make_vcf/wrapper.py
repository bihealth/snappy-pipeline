# -*- coding: utf-8 -*-
"""Wrapper for running Melt make_vcf
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

genotype_dir=$(dirname {snakemake.input.genotype} | head -n 1)
ls $genotype_dir/*.{snakemake.wildcards.me_type}.tsv \
| sort \
> {snakemake.output.list_txt}

# TODO: allowing 100% no-call wise?
melt_mei -Xmx2G MakeVCF \
    -genotypingdir $genotype_dir \
    -h {snakemake.config[static_data_config][reference][path]} \
    -j 100 \
    -t $(melt_mei_path)/me_refs/{snakemake.config[step_config][wgs_mei_calling][melt][me_refs_infix]}/{snakemake.wildcards.me_type}_MELT.zip \
    -p $(dirname {snakemake.input.group_analysis}) \
    -w $(dirname {snakemake.output.done}) \
    -o $(dirname {snakemake.output.done})

# Make file more VCF conforming
perl -p -e 's/ID=GL,Number=3/ID=GL,Number=G/' $(dirname {snakemake.output.done})/{snakemake.wildcards.me_type}.final_comp.vcf \
| bgzip -c \
> {snakemake.output.vcf}
tabix -f {snakemake.output.vcf}
"""
)
