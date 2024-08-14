# -*- coding: utf-8 -*-
"""Wrapper for running bcftools filter"""

from snakemake.shell import shell
from snakemake.script import snakemake

args = snakemake.params["args"]
filter_name = args["filter_name"]
expression = (
    '--include "{}"'.format(args["include"])
    if "include" in args
    else '--exclude "{}"'.format(args["exclude"])
)

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

bcftools filter --soft-filter {filter_name} --mode + \
    {expression} \
    -O z -o {snakemake.output.vcf} \
    {snakemake.input.vcf}
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
