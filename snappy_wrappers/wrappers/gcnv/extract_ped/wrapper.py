# -*- coding: utf-8 -*-
# Extract pedigree members only from VCF file.

from snakemake.shell import shell

shell(
    r"""
# -----------------------------------------------------------------------------
# Redirect stderr to log file by default and enable printing executed commands
exec &> >(tee -a "{snakemake.log}")
set -x
# -----------------------------------------------------------------------------

set -euo pipefail

export TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT

echo '{snakemake.params.ped_members}' \
| tr ' ' '\n' \
| LANG=C sort \
>$TMPDIR/samples.txt

if grep . $TMPDIR/samples.txt; then
    bcftools view \
        --force-samples \
        --samples-file $TMPDIR/samples.txt \
        --output-type u \
        {snakemake.input.vcf} \
    | bcftools view \
        --output-file {snakemake.output.vcf} \
        --output-type z \
        --include '(GT == "alt")'
else
    # no sample, just keep header
    bcftools view \
        --force-samples \
        --samples-file $TMPDIR/samples.txt \
        {snakemake.input.vcf} \
    | zgrep '^#' \
    | bgzip -c \
    > {snakemake.output.vcf}
fi

tabix -f {snakemake.output.vcf}

pushd $(dirname {snakemake.output.vcf})
md5sum $(basename {snakemake.output.vcf}) >$(basename {snakemake.output.vcf}).md5
md5sum $(basename {snakemake.output.vcf}).tbi >$(basename {snakemake.output.vcf}).tbi.md5
"""
)
