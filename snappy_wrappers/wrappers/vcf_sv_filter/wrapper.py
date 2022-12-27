# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for vcf_sv_filter

- Remove LowQual variants
- Add various annotations
"""

import os

from snakemake import shell

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

shell.executable("/bin/bash")

base_dir = os.path.dirname(os.path.realpath(__file__))

shell(
    r"""
set -x

export TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT

# Load library with helper functions.
source {base_dir}/../wgs_sv_filtration/funcs.sh

# Run through VCF SV filter and limit to samples in family.

set +e
samples=$(
    samples_vcf_ped {snakemake.input.sv_bcf} {snakemake.input.ped} \
    | tr '\n' ',' \
    | sed -e 's/,$//'
)
set -e

if [[ {snakemake.wildcards.caller} == delly2 ]]; then
    bcftools view --threads 4 --force-samples -s $samples -e 'FILTER == "LowQual"' {snakemake.input.sv_bcf}
else
    bcftools view --threads 4 --force-samples -s $samples                          {snakemake.input.sv_bcf}
fi \
| bcftools query -f "%ID\n" -i 'GT ~ "1"' \
> $TMPDIR/ids.txt

( \
    bcftools view -h {snakemake.input.sv_bcf};
    bcftools view --threads 4 -H {snakemake.input.sv_bcf} \
    | grep -w -F -f $TMPDIR/ids.txt \
) \
| bgzip --threads 4 -c \
> $TMPDIR/tmp.vcf.gz

tabix -f $TMPDIR/tmp.vcf.gz

# Run through VCF SV filter

time python3 {base_dir}/vcf_sv_filter.py \
    --ped-file {snakemake.input.ped} \
    --input-vcf $TMPDIR/tmp.vcf.gz \
    --small-var-vcf {snakemake.input.var_vcf} \
    $(if [[ -n "{snakemake.config[step_config][wgs_sv_annotation][path_alu_bed]}" ]]; then echo --alu-bed {snakemake.config[step_config][wgs_sv_annotation][path_alu_bed]}; fi) \
    $(if [[ -n "{snakemake.config[step_config][wgs_sv_annotation][path_db_bed]}" ]]; then echo --db-bed {snakemake.config[step_config][wgs_sv_annotation][path_db_bed]}; fi) \
| bgzip -c \
> $TMPDIR/tmp2.vcf.gz

bcftools view \
    --threads 4 \
    --force-samples \
    -s $samples \
    $TMPDIR/tmp2.vcf.gz \
| bcftools view \
    -i 'GT ~ "1"' \
    -O z \
    -o {snakemake.output.vcf}

tabix -f {snakemake.output.vcf}

pushd $(dirname {snakemake.output.vcf})
md5sum $(basename {snakemake.output.vcf}) >$(basename {snakemake.output.vcf}).md5
md5sum $(basename {snakemake.output.vcf_tbi}) >$(basename {snakemake.output.vcf_tbi}).md5
"""
)
