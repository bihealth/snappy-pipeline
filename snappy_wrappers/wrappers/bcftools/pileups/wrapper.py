# -*- coding: utf-8 -*-
"""Wrapper for running bcftools mpileup"""

from typing import TYPE_CHECKING

from snakemake.shell import shell

if TYPE_CHECKING:
    from snakemake.script import snakemake


args = snakemake.params["args"]
reference_path = args["reference_path"]
max_depth = args["max_depth"]

# FIXME: "locii" only ever gets set as the input, never as a parameter in args
if intervals := args["intervals"]:
    locii = "-r " + intervals
elif "locii" in snakemake.input.keys():
    locii = "-R " + snakemake.input.locii
elif locii_arg := args.get("locii"):
    locii = "-R " + locii_arg
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
    --max-depth {max_depth} \
    -f {reference_path} \
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
