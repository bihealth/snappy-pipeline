# -*- coding: utf-8 -*-
"""Wrapper for finding heterozygous variants with bcftools"""

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

# Convert minimum B-allele fraction into ratio of alternative to reference alleles
min_ratio = config["min_baf"] / (1 - config["min_baf"])
max_ratio = 1 / min_ratio

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

only_one_variant="N_ALT=2 & FORMAT/AD[:2]=0"
min_depth="FORMAT/AD[:0]>{config[min_depth]} & FORMAT/AD[:1]>{config[min_depth]}"
hetero="{min_ratio}*FORMAT/AD[:0]<=FORMAT/AD[:1] & FORMAT/AD[:1]<={max_ratio}*FORMAT/AD[:0]"

bcftools mpileup \
    {locii} \
    --max-depth {config[max_depth]} \
    -f {snakemake.config[static_data_config][reference][path]} \
    -a "FORMAT/AD" \
    {snakemake.input.bam} \
    | bcftools filter \
    --include "$only_one_variant & $min_depth & $hetero" \
    -O z -o {snakemake.output.vcf}
tabix {snakemake.output.vcf}

pushd $(dirname {snakemake.output.vcf})
md5sum $(basename {snakemake.output.vcf}) > $(basename {snakemake.output.vcf_md5})
md5sum $(basename {snakemake.output.vcf_tbi}) > $(basename {snakemake.output.vcf_tbi_md5})
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
