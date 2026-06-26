# -*- coding: utf-8 -*-
"""
Wrapper for quantify gene and transcript expression
"""

from snakemake.shell import shell

__author__ = "Pham Gia Cuong"
__email__ = "pham.gia-cuong@bih-charite.de"

step = snakemake.config["pipeline_step"]["name"]
config = snakemake.config["step_config"][step]

shell(
    r"""
# -----------------------------------------------------------------------------
# Redirect stderr to log file by default and enable printing executed commands
exec 2> >(tee -a "{snakemake.log.log}")
set -x
# -----------------------------------------------------------------------------

# Write out information about conda installation
conda list > {snakemake.log.conda_list}
conda info > {snakemake.log.conda_info}


stringtie \
    -G {snakemake.config[static_data_config][features][path]} \
    -p 4 \
    -v \
    -o {snakemake.output.expression} \
    {snakemake.input.rna_bam}

pushd $(dirname {snakemake.output.expression})
md5sum $(basename {snakemake.output.expression}) > $(basename {snakemake.output.expression_md5})
"""
)


# Compute MD5 sums of logs
shell(
    r"""
md5sum {snakemake.log.log} > {snakemake.log.log_md5}
md5sum {snakemake.log.conda_list} > {snakemake.log.conda_list_md5}
md5sum {snakemake.log.conda_info} > {snakemake.log.conda_info_md5}
"""
)
