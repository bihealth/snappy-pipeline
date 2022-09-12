# -*- coding: utf-8 -*-
"""Wrapper for running Delly2's CNV re-genotyping step"""

from snakemake.shell import shell

shell(
    r"""
# -----------------------------------------------------------------------------
# Redirect stderr to log file by default and enable printing executed commands
exec &> >(tee -a "{snakemake.log.log}")
set -x
# -----------------------------------------------------------------------------

# Write out information about conda installation
conda list > {snakemake.log.conda_list}
conda info > {snakemake.log.conda_info}
md5sum {snakemake.log.conda_list} > {snakemake.log.conda_list_md5}
md5sum {snakemake.log.conda_info} > {snakemake.log.conda_info_md5}


delly cnv --segmentation \
    --mappability {snakemake.config[step_config][wgs_cnv_calling][delly2][mappability]} \
    --genome {snakemake.config[static_data_config][reference][path]} \
    --vcffile {snakemake.input.bcf} \
    --outfile {snakemake.output.bcf} \
    {snakemake.input.bam}

tabix -f {snakemake.output.bcf}

pushd $(dirname {snakemake.output.bcf})
md5sum $(basename {snakemake.output.bcf}) > $(basename {snakemake.output.bcf_md5})
md5sum $(basename {snakemake.output.csi}) > $(basename {snakemake.output.csi_md5})
"""
)

# Compute MD5 sums of logs
shell(
    r"""
md5sum {snakemake.log.log} > {snakemake.log.log_md5}
"""
)
