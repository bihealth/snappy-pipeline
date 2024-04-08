# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper for filtering WGS SV results for mode of inheritance."""

# TODO: works for trios, singletons but NOT FOR MORE COMPLICATED CASES
# TODO: how to make work if only one parent present?

import os

from snakemake import shell

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

# shell.executable('/bin/bash') # XXX

base_dir = os.path.dirname(os.path.realpath(__file__))

shell(
    r"""
set -x

export TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT

# Load library with helper functions.
source {base_dir}/../funcs.sh

# Get name and number of index, father, and mother ------------------------------------------------

index={snakemake.wildcards.index_library}
father=$(awk '($2 == "R13_1-N1-DNA1-WGS1") {{ print $3; }}' {snakemake.input.ped})
mother=$(awk '($2 == "R13_1-N1-DNA1-WGS1") {{ print $4; }}' {snakemake.input.ped})

index_no=$(get_index {snakemake.input.vcf} "$index")
father_no=$(get_index {snakemake.input.vcf} "$father")
mother_no=$(get_index {snakemake.input.vcf} "$mother")

# Definition of the filtration functions ----------------------------------------------------------

# Dominant/de novo/recessive hom can be solved with a simple filter expression.
simple_filter() {{
    bcftools view \
        -i "$1" \
        -O z \
        -o {snakemake.output.vcf} \
        {snakemake.input.vcf}
}}

# Recessive het. comp. in the case of only one parent
rec_hc_one_parent() {{
    # Number of parents present in pedigree
    parents=0
    [[ ! -z "$1" ]] && parents+=1
    [[ ! -z "$2" ]] && parents+=1

    # One or both parents missing, give up on being smart.
    bcftools view \
        -i "(UNAFFECTED_CARRIERS == $parents) && (BACKGROUND_CARRIERS == 0)" \
        -O z \
        -o {snakemake.output.vcf} \
        {snakemake.input.vcf}
}}

# Perform filtration for ARHC variants (2SVs affecting the same tad, one from father and one
# from mother)
rec_hc_two_parents() {{
    father_no=$1
    mother_no=$2

    tads_bed={snakemake.config[step_config][wgs_sv_filtration][region_beds][all_tads]}

    # Create mother SV file
    exp_mother="(UNAFFECTED_CARRIERS == 1)"
    exp_mother+=" && (BACKGROUND_CARRIERS == 0)"
    exp_mother+=" && (GT[$mother_no] != \"alt\")"

    bcftools view \
        -i "$exp_mother" \
        -O z \
        -o $TMPDIR/mother_SVs.vcf.gz \
        {snakemake.input.vcf}

    # Create father SV file
    exp_father="(UNAFFECTED_CARRIERS == 1)"
    exp_father+=" && (BACKGROUND_CARRIERS == 0)"
    exp_father+=" && (GT[$father_no] == \"alt\")"

    bcftools view \
        -i "$exp_father" \
        -O z \
        -o $TMPDIR/father_SVs.vcf.gz \
        {snakemake.input.vcf}

    # Perform interval intersections to derive final AR het.comp. set of SVs.
    bedtools intersect -wa -a $tads_bed -b $TMPDIR/mother_SVs.vcf.gz \
    | sort -k1,1V -k2,2n \
    | uniq \
    > $TMPDIR/mother_SVs.tads.bed

    bedtools intersect -wa -a $tads_bed -b $TMPDIR/father_SVs.vcf.gz \
    | sort -k1,1V -k2,2n \
    | uniq \
    > $TMPDIR/father_SVs.tads.bed

    sort $TMPDIR/mother_SVs.tads.bed \
    > $TMPDIR/mother_SVs.tads.sorted.bed

    sort $TMPDIR/father_SVs.tads.bed \
    > $TMPDIR/father_SVs.tads.sorted.bed

    comm -12 $TMPDIR/mother_SVs.tads.sorted.bed $TMPDIR/father_SVs.tads.sorted.bed \
    > $TMPDIR/commom_tads.bed

    sort -k1,1 -k2,2n $TMPDIR/commom_tads.bed \
    > $TMPDIR/commom_tads.sorted.bed

    bedtools intersect -header -wa -a $TMPDIR/mother_SVs.vcf.gz -b $TMPDIR/commom_tads.sorted.bed \
    | sort -k1,1V -k2,2n \
    | uniq \
    > $TMPDIR/ARHC.frommother.vcf

    bedtools intersect -header -wa -a $TMPDIR/father_SVs.vcf.gz -b $TMPDIR/commom_tads.sorted.bed \
    | sort -k1,1V -k2,2n \
    | uniq \
    > $TMPDIR/ARHC.fromfather.vcf

    bgzip $TMPDIR/ARHC.frommother.vcf
    bgzip $TMPDIR/ARHC.fromfather.vcf

    tabix -f $TMPDIR/ARHC.frommother.vcf.gz
    tabix -f $TMPDIR/ARHC.fromfather.vcf.gz

    bcftools concat -a -O z -o {snakemake.output.vcf} \
    $TMPDIR/ARHC.frommother.vcf.gz $TMPDIR/ARHC.fromfather.vcf.gz
}}

# Actual filtration -------------------------------------------------------------------------------

## TODO: think whether to generate two files (father / mother) adapt the snakemake.wilcards.inheritance somewhere
## or one and filter for father/mother at the ARHC step

case "{snakemake.wildcards.inheritance}" in
    dominant)
        # Perform filtration for variants dominantly transmitted by father or mother
        simple_filter '(UNAFFECTED_CARRIERS == 1) && (BACKGROUND_CARRIERS == 0)'
        ;;
    de_novo)
        # Perform filtration for de novo variants
        simple_filter '(UNAFFECTED_CARRIERS == 0) && (BACKGROUND_CARRIERS == 0)'
        ;;
    recessive_hom)
        # Perform filtration for homozygous variant present in both parents
        simple_filter '(UNAFFECTED_CARRIERS == 2) && (BACKGROUND_CARRIERS == 0)'
        ;;
    recessive_hc)
        if [[ -z "$father_no" ]] || [[ -z "$mother_no" ]]; then
            rec_hc_one_parent "$father_no" "$mother_no"
        else
            rec_hc_two_parents "$father_no" "$mother_no"
        fi
    ;;
    *)  # else, "all"
    cp {snakemake.input.vcf} {snakemake.output.vcf}
esac

tabix -f {snakemake.output.vcf}

# Compute checksums -------------------------------------------------------------------------------

pushd $(dirname {snakemake.output.vcf})
md5sum $(basename {snakemake.output.vcf}) >$(basename {snakemake.output.vcf}).md5
md5sum $(basename {snakemake.output.vcf_tbi}) >$(basename {snakemake.output.vcf_tbi}).md5
"""
)
