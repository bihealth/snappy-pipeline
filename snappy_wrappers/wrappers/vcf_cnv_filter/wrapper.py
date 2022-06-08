# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for soft-annotating CNV VCF files

- Add annotations based on background.
"""

import os

from snakemake import shell

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

shell.executable("/bin/bash")

# Get path to this file's (wrapper.py) directory.
base_dir = os.path.dirname(os.path.realpath(__file__))

bed_files = []
for bed_file in snakemake.config["step_config"][snakemake.params.step_name].get("bed_files", []):
    bed_files.append(
        "{}|{}|{}".format(bed_file["name"], bed_file["description"], bed_file["path"]).replace(
            " ", "_"
        )
    )
bed_files = " ".join(map(lambda x: repr(x), bed_files))

# TODO: implement support for more than trios

shell(
    r"""
set -x

## Write out information about conda installation.
#conda list >{{snakemake.log.conda_list}}
#conda info >{{snakemake.log.conda_info}}
#md5sum {{snakemake.log.conda_list}} >{{snakemake.log.conda_list_md5}}
#md5sum {{snakemake.log.conda_info}} >{{snakemake.log.conda_info_md5}}
#
## Also pipe stderr to log file
#if [[ -n "{{snakemake.log.log}}" ]]; then
#    if [[ "$(set +e; tty; set -e)" != "" ]]; then
#        rm -f "{{snakemake.log.log}}" && mkdir -p $(dirname {{snakemake.log.log}})
#        exec 2> >(tee -a "{{snakemake.log.log}}" >&2)
#    else
#        rm -f "{{snakemake.log.log}}" && mkdir -p $(dirname {{snakemake.log.log}})
#        echo "No tty, logging disabled" >"{{snakemake.log.log}}"
#    fi
#fi

export TMPDIR=$(mktemp -d)
# trap "rm -rf $TMPDIR" EXIT

# Load library with helper functions.
source {base_dir}/../wgs_sv_filtration/funcs.sh

# Get name and number of index, father, and mother.
index={snakemake.wildcards.index_ngs_library}
father=$(awk '($2 == "'$index'") {{ print $3; }}' {snakemake.input.ped})
mother=$(awk '($2 == "'$index'") {{ print $4; }}' {snakemake.input.ped})

index_no=$(get_index {snakemake.input.vcf} "$index")
father_no=$(get_index {snakemake.input.vcf} "$father")
mother_no=$(get_index {snakemake.input.vcf} "$mother")

include="(GT[$index_no] == \"alt\")"
test -n "$father_no" && include+="|| (GT[$father_no] == \"alt\")"
test -n "$mother_no" && include+="|| (GT[$mother_no] == \"alt\")"

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

set +e
samples=$(
    samples_vcf_ped {snakemake.input.vcf} {snakemake.input.ped} \
    | tr '\n' ',' \
    | sed -e 's/,$//'
)
set -e

if [[ -n "$samples" ]]; then
    time python3 {base_dir}/vcf_cnv_filter.py \
        --ped-file {snakemake.input.ped} \
        --input-vcf $TMPDIR/tmp.vcf.gz \
        $(for bed_file in {bed_files}; do \
            echo --annotation-bed $bed_file
        done) \
    | bcftools view \
        -s "$samples" \
        -O z \
        -o {snakemake.output.vcf}
else
    bcftools view -s "$samples" --force-samples $TMPDIR/tmp.vcf.gz \
    | grep '^#' \
    > $TMPDIR/tmp2.vcf

    echo '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">' \
    > $TMPDIR/header.txt

    bcftools annotate \
        -h $TMPDIR/header.txt \
        -Oz \
        -o {snakemake.output.vcf} \
        $TMPDIR/tmp2.vcf
fi

tabix -f {snakemake.output.vcf}

pushd $(dirname {snakemake.output.vcf})
md5sum $(basename {snakemake.output.vcf}) >$(basename {snakemake.output.vcf}).md5
md5sum $(basename {snakemake.output.tbi}) >$(basename {snakemake.output.tbi}).md5
"""
)
