# -*- coding: utf-8 -*-
"""Wrapper for running bcftools mpileup
"""

from snakemake.shell import shell

step = snakemake.config["pipeline_step"]["name"]
config = snakemake.config["step_config"][step]

if "args" in snakemake.params and "intervals" in snakemake.params["args"]:
    locii = "-r " + snakemake.params["args"]["intervals"]
elif "locii" in snakemake.input.keys():
    locii = "-R " + snakemake.input.locii
elif "locii" in config and config["locii"]:
    locii = "-R " + config["locii"]
else:
    locii = ""

# Actually run the script.
shell(
    r"""
# -----------------------------------------------------------------------------
# Redirect stderr to log file by default and enable printing executed commands
exec &> >(tee -a "{snakemake.log.log}")
set -x
# -----------------------------------------------------------------------------
export TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT

# Write out information about conda installation
conda list > {snakemake.log.conda_list}
conda info > {snakemake.log.conda_info}

bcftools mpileup \
    {locii} \
    --max-depth {config[max_depth]} \
    -f {snakemake.config[static_data_config][reference][path]} \
    -a "FORMAT/AD" \
    -O z -o {snakemake.output.vcf} \
    {snakemake.input.bam}
tabix {snakemake.output.vcf}

pushd $(dirname {snakemake.output.vcf})
md5sum $(basename {snakemake.output.vcf}) > $(basename {snakemake.output.vcf_md5})
md5sum $(basename {snakemake.output.vcf_tbi}) > $(basename {snakemake.output.vcf_tbi_md5})
popd
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
