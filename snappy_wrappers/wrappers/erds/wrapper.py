# -*- coding: utf-8 -*-
"""Wrapper for running ERDS call step
"""

from snakemake.shell import shell

__author__ = "Manuel Holtgrewe"
__email__ = "manuel.holtgrewe@bih-charite.de"

shell(
    r"""
# -----------------------------------------------------------------------------
# Redirect stderr to log file by default and enable printing executed commands
exec &> >(tee -a "{snakemake.log}")
set -x
# -----------------------------------------------------------------------------

export TMPDIR=$(mktemp -d)
# TODO: cleanup after us!
trap "rm -rf $TMPDIR" EXIT

mkdir $TMPDIR/{{out,tmp}}

# Filter VCF file to genotypes occuring in the sample at hand.
bcftools view \
    -O u \
    -s "{snakemake.wildcards.library_name}" \
    {snakemake.input.vcf} \
| bcftools view \
    -i 'GT ~ "1"' \
    -O v \
    -o $TMPDIR/tmp/tmp.vcf

# Run ERDS Pipeline.
erds_pipeline \
    -b {snakemake.input.bam} \
    -v $TMPDIR/tmp/tmp.vcf \
    -o $TMPDIR/out \
    -r {snakemake.config[static_data_config][reference][path]}

# # Completely read in VCF file to check for errors.
# set -e
# python3 {{base_dir}}/check_vcf.py \
#     "$TMPDIR/out/{snakemake.wildcards.library_name}.erds.vcf"

# Convert result to output file, and create tabix index.
bgzip -c "$TMPDIR/out/{snakemake.wildcards.library_name}.erds.vcf" \
> {snakemake.output.vcf}
tabix -f {snakemake.output.vcf}

# Finally, create MD5 sum files.
pushd $(dirname {snakemake.output.vcf})
md5sum $(basename {snakemake.output.vcf}) >$(basename {snakemake.output.vcf}).md5
md5sum $(basename {snakemake.output.vcf}).tbi >$(basename {snakemake.output.vcf}).tbi.md5
"""
)
