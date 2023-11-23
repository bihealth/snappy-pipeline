# -*- coding: utf-8 -*-
"""Wrapper for running bcftools mpileup
"""

from snakemake.shell import shell

step = snakemake.config["pipeline_step"]["name"]
config = snakemake.config["step_config"][step]

if "include" in config or "exclude" in config:
    include = config.get("include", "")
    exclude = config.get("exclude", "")
elif "bcftools" in config and ("include" in config["bcftools"] or "exclude" in config["bcftools"]):
    include = config["bcftools"].get("include", "")
    exclude = config["bcftools"].get("exclude", "")
elif "filter_nb" in snakemake.wildcards.keys():
    filter_nb = int(snakemake.wildcards["filter_nb"]) - 1
    include = config["filter_list"][filter_nb]["bcftools"].get("include", "")
    exclude = config["filter_list"][filter_nb]["bcftools"].get("exclude", "")
else:
    include = ""
    exclude = ""
if include:
    include = '--include "' + include + '"'
if exclude:
    exclude = '--exclude "' + exclude + '"'

if "filter_name" in config:
    filter_name = config.get("filter_name", "")
elif "bcftools" in config and "filter_name" in config["bcftools"]:
    filter_name = config["bcftools"].get("filter_name", "")
elif "filter_nb" in snakemake.wildcards.keys():
    filter_name = "bcftools_{}".format(int(snakemake.wildcards["filter_nb"]))
else:
    filter_name = "+"

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
    {include} {exclude} \
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
