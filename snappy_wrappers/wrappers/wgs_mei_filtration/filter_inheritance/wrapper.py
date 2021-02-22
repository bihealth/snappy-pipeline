# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for inheritance filter for wgs_mei_filtration.
"""

# TODO: works for trios, singletons but NOT FOR MORE COMPLICATED CASES
# TODO: how to make work if only one parent present?
# TODO: consolidate with same filter for SV

import os

from snakemake import shell

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>"

shell.executable("/bin/bash")

base_dir = os.path.dirname(os.path.realpath(__file__))

shell(
    r"""
set -x

export TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT

# Load library with helper functions.
source {base_dir}/../../wgs_sv_filtration/funcs.sh

# Get name and number of index, father, and mother.
index={snakemake.wildcards.index_library}
father=$(awk '($2 == "R13_1-N1-DNA1-WGS1") {{ print $3; }}' {snakemake.input.ped})
mother=$(awk '($2 == "R13_1-N1-DNA1-WGS1") {{ print $4; }}' {snakemake.input.ped})

index_no=$(get_index {snakemake.input.vcf} "$index")
father_no=$(get_index {snakemake.input.vcf} "$father")
mother_no=$(get_index {snakemake.input.vcf} "$mother")

# Perform the actual filtration
tads_bed={snakemake.config[step_config][wgs_mei_filtration][region_beds][all_tads]}

## TODO: think whether to generate two files (father / mother) adapt the snakemake.wilcards.inheritance somewhere
## or one and filter for father/mother at the ARHC step

if [[ "{snakemake.wildcards.inheritance}" == dominant ]]; then
    # Perform filtration for variants dominantly transmitted by father or mother
    bcftools view \
        -i '(UNAFFECTED_CARRIERS == 1) && (BACKGROUND_CARRIERS == 0)' \
        -O z \
        -o {snakemake.output.vcf} \
        {snakemake.input.vcf}
elif [[ "{snakemake.wildcards.inheritance}" == de_novo ]]; then
    # Perform filtration for de novo variants
    bcftools view \
        -i '(UNAFFECTED_CARRIERS == 0) && (BACKGROUND_CARRIERS == 0)' \
        -O z \
        -o {snakemake.output.vcf} \
        {snakemake.input.vcf}
elif [[ "{snakemake.wildcards.inheritance}" == recessive_hom ]]; then
    # Perform filtration for homozygous variant present in both parents
    bcftools view \
        -i '(UNAFFECTED_CARRIERS == 2) && (BACKGROUND_CARRIERS == 0)' \
        -O z \
        -o {snakemake.output.vcf} \
        {snakemake.input.vcf}
elif [[ "{snakemake.wildcards.inheritance}" == recessive_hc ]]; then
    if [[ -z "$father_no" ]] || [[ -z "$father_no" ]]; then
        # Number of parents present in pedigree
        parents=0
        [[ ! -z "$father_no" ]] && parents+=1
        [[ ! -z "$mother_no" ]] && parents+=1

        # One or both parents missing, give up on being smart.
        bcftools view \
            -i "(UNAFFECTED_CARRIERS == $parents) && (BACKGROUND_CARRIERS == 0)" \
            -O z \
            -o {snakemake.output.vcf} \
            {snakemake.input.vcf}
    else
        # Perform filtration for ARHC variants (2MEIs affecting the same tad, one from father and one
        # from mother)

        # Create mother MEI file
        exp_mother="(UNAFFECTED_CARRIERS == 1)"
        exp_mother+=" && (BACKGROUND_CARRIERS == 0)"
        exp_mother+=" && (GT[$mother_no] != \"0/0\")"

        bcftools view \
            -i "$exp_mother" \
            -O z \
            -o $TMPDIR/mother_meis.vcf.gz \
            {snakemake.input.vcf}

        # Create father MEI file
        exp_father="(UNAFFECTED_CARRIERS == 1)"
        exp_father+=" && (BACKGROUND_CARRIERS == 0)"
        exp_father+=" && (GT[$father_no] != \"0/0\")"

        bcftools view \
            -i "$exp_father" \
            -O z \
            -o $TMPDIR/father_meis.vcf.gz \
            {snakemake.input.vcf}

        # Perform interval intersections to derive final AR het.comp. set of MEIs.
        bedtools intersect -wa -a $tads_bed -b $TMPDIR/mother_meis.vcf.gz \
        | sort -k1,1V -k2,2n \
        | uniq \
        > $TMPDIR/mother_meis.tads.bed

        bedtools intersect -wa -a $tads_bed -b $TMPDIR/father_meis.vcf.gz \
        | sort -k1,1V -k2,2n \
        | uniq \
        > $TMPDIR/father_meis.tads.bed

        sort $TMPDIR/mother_meis.tads.bed \
        > $TMPDIR/mother_meis.tads.sorted.bed

        sort $TMPDIR/father_meis.tads.bed \
        > $TMPDIR/father_meis.tads.sorted.bed

        comm -12 $TMPDIR/mother_meis.tads.sorted.bed $TMPDIR/father_meis.tads.sorted.bed \
        > $TMPDIR/commom_tads.bed

        sort -k1,1 -k2,2n $TMPDIR/commom_tads.bed \
        > $TMPDIR/commom_tads.sorted.bed

        bedtools intersect -header -wa -a $TMPDIR/mother_meis.vcf.gz -b $TMPDIR/commom_tads.sorted.bed \
        | sort -k1,1V -k2,2n \
        | uniq \
        > $TMPDIR/ARHC.frommother.vcf

        bedtools intersect -wa -a $TMPDIR/father_meis.vcf.gz -b $TMPDIR/commom_tads.sorted.bed \
        | sort -k1,1V -k2,2n \
        | uniq \
        > $TMPDIR/ARHC.fromfather.vcf

        grep ^\# $TMPDIR/ARHC.frommother.vcf \
        > $TMPDIR/ARHC.vcf

        cat $TMPDIR/ARHC.frommother.vcf $TMPDIR/ARHC.fromfather.vcf \
        | {{ grep -v ^\# || true; }} \
        | sort -k1,1V -k2,2n \
        >> $TMPDIR/ARHC.vcf

        cat $TMPDIR/ARHC.vcf \
        | bgzip -c \
        > $TMPDIR/ARHC.vcf.gz
        cp $TMPDIR/ARHC.vcf.gz {snakemake.output.vcf}
    fi
else  # else, "all"
    cp {snakemake.input.vcf} {snakemake.output.vcf}
fi

tabix -f {snakemake.output.vcf}

pushd $(dirname {snakemake.output.vcf})
md5sum $(basename {snakemake.output.vcf}) >$(basename {snakemake.output.vcf}).md5
md5sum $(basename {snakemake.output.tbi}) >$(basename {snakemake.output.tbi}).md5
"""
)
