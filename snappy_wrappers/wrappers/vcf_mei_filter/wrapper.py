# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for soft-annotating MEI VCF files

- Add annotations based on background.
"""

import os

from snakemake import shell

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

shell.executable("/bin/bash")

base_dir = os.path.dirname(os.path.realpath(__file__))

# Get path to this file's (wrapper.py) directory.
base_dir = os.path.dirname(os.path.realpath(__file__))

# TODO: implement support for more than trios

shell(
    r"""
set -x

export TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT

# Load library with helper functions.
source {base_dir}/../wgs_sv_filtration/funcs.sh

# Get name and number of index, father, and mother.
index={snakemake.wildcards.index_ngs_library}
father=$(awk '($2 == "'$index'") {{ print $3; }}' {snakemake.input.ped})
mother=$(awk '($2 == "'$index'") {{ print $4; }}' {snakemake.input.ped})

index_no=$(get_index {snakemake.input.vcf} "$index")
father_no=$(get_index {snakemake.input.vcf} "$father")
mother_no=$(get_index {snakemake.input.vcf} "$mother")

include="(GT[$index_no] ~ \"1\")"
test -n "$father_no" && include+="|| (GT[$father_no] ~ \"1\")"
test -n "$mother_no" && include+="|| (GT[$mother_no] ~ \"1\")"

# Limit to variants present in family.
bcftools view \
    --threads 4 \
    --force-samples \
    -i "$include" \
    {snakemake.input.vcf} \
| bgzip -c \
> $TMPDIR/tmp.vcf.gz

tabix -f $TMPDIR/tmp.vcf.gz

# Run through VCF SV filter and limit to samples in family.

time python3 {base_dir}/vcf_mei_filter.py \
    --ped-file {snakemake.input.ped} \
    --input-vcf $TMPDIR/tmp.vcf.gz \
| bcftools view \
    -s "$(
            samples_vcf_ped {snakemake.input.vcf} {snakemake.input.ped} \
            | tr '\n' ',' \
            | sed -e 's/,$//'
        )" \
    -O z \
    -o {snakemake.output.vcf}

tabix -f {snakemake.output.vcf}

pushd $(dirname {snakemake.output.vcf})
md5sum $(basename {snakemake.output.vcf}) >$(basename {snakemake.output.vcf}).md5
md5sum $(basename {snakemake.output.tbi}) >$(basename {snakemake.output.tbi}).md5
"""
)
