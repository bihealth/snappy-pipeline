#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

snappy-add_ped_header_usage() {
    >&2 echo Usage: add_ped_header PED VCF OUT.VCF
    exit 1
}

snappy-add_ped_header()
{
    PED=$1
    VCF=$2
    OUT=$3

    VCF_PED=$(mktemp)
    VCF_HEADER=$(mktemp)

    test -z $PED && snappy-add_ped_header_usage
    test -z $VCF && snappy-add_ped_header_usage
    test -z $OUT && snappy-add_ped_header_usage

    python $DIR/ped_to_vcf_header.py --ped-file $PED --output $VCF_PED
    bcftools view --header-only $VCF > $VCF_HEADER
    bcftools reheader \
        --header <(sed "/#CHROM/e cat $VCF_PED" $VCF_HEADER) \
        --output $OUT \
        $VCF
}
