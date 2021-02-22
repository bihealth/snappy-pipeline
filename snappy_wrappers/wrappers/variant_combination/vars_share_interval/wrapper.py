# -*- coding: utf-8 -*-
"""Snakemake wrapper for combination of variants.
"""

import os

from snakemake import shell

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>"

shell.executable("/bin/bash")

step_config = snakemake.config["step_config"]["variant_combination"]
combinations = {}
for combination in step_config["combinations"]:
    combinations[combination["name"]] = combination

args = snakemake.params["args"]
intervals_bed = args["intervals_bed"]
left_source = combinations[snakemake.wildcards.combination]["left"].split(":", 1)[0]
right_source = combinations[snakemake.wildcards.combination]["right"].split(":", 1)[0]

# Get path to this file's (wrapper.py) directory.
base_dir = os.path.dirname(os.path.realpath(__file__))

shell(
    r"""
>&2 echo INTERVALS: {intervals_bed}
>&2 echo LEFT SOURCE: {left_source}
>&2 echo RIGHT SOURCE: {right_source}

set -x
set -euo pipefail

# Setup auto-cleaned TMPDIR
export TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT

# Load library with helper functions.
source {base_dir}/../../wgs_sv_filtration/funcs.sh

# Get name and number of index, father, and mother.
index={snakemake.wildcards.index_library}
father=$(awk '($2 == "'$index'") {{ print $3; }}' {snakemake.input.ped})
mother=$(awk '($2 == "'$index'") {{ print $4; }}' {snakemake.input.ped})

left_index_no=$(get_index {snakemake.input.left_vcf} "$index")
left_father_no=$(get_index {snakemake.input.left_vcf} "$father")
left_mother_no=$(get_index {snakemake.input.left_vcf} "$mother")

right_index_no=$(get_index {snakemake.input.right_vcf} "$index")
right_father_no=$(get_index {snakemake.input.right_vcf} "$father")
right_mother_no=$(get_index {snakemake.input.right_vcf} "$mother")

# Short-circuit in the case of missing father/mother
if [[ -z "$left_father_no" ]] || [[ -z "$left_mother_no" ]] || \
        [[ -z "$right_father_no" ]] || [[ -z "$right_mother_no" ]]; then
    touch \
        {snakemake.output.vcf} \
        {snakemake.output.vcf}.tbi
    pushd $(dirname {snakemake.output.vcf})
    md5sum $(basename {snakemake.output.vcf}) >$(basename {snakemake.output.vcf}).md5
    md5sum $(basename {snakemake.output.tbi}) >$(basename {snakemake.output.tbi}).md5
    exit 0
fi

# We will create temporary files that contain only the index, mother and father data,
# in this order.
bcftools view -O z -o $TMPDIR/SVs.unfiltered.vcf.gz -s "$index,$mother,$father" \
    {snakemake.input.left_vcf}
bcftools view -O z -o $TMPDIR/SNVs.unfiltered.vcf.gz -s "$index,$mother,$father" \
    {snakemake.input.right_vcf}

# Create mother files
bcftools view \
    -i "GT[1] != \"0/0\"" \
    -O z \
    -o $TMPDIR/SVs.mother.vcf.gz \
    $TMPDIR/SVs.unfiltered.vcf.gz
bcftools view \
    -i "GT[1] != \"0/0\"" \
    -O z \
    -o $TMPDIR/SNVs.mother.vcf.gz \
    $TMPDIR/SNVs.unfiltered.vcf.gz

# NB: does nothing for now because these are the hom files
echo $(zcat $TMPDIR/SVs.unfiltered.vcf.gz | grep -v ^\# | wc -l)
echo $(zcat $TMPDIR/SVs.mother.vcf.gz | grep -v ^\# | wc -l)
echo $(zcat $TMPDIR/SNVs.unfiltered.vcf.gz | grep -v ^\# | wc -l)
echo $(zcat $TMPDIR/SNVs.mother.vcf.gz | grep -v ^\# | wc -l)

# Create father files
bcftools view \
    -i "GT[2] != \"0/0\"" \
    -O z \
    -o $TMPDIR/SVs.father.vcf.gz \
    $TMPDIR/SVs.unfiltered.vcf.gz
bcftools view \
    -i "GT[2] != \"0/0\"" \
    -O z \
    -o $TMPDIR/SNVs.father.vcf.gz \
    $TMPDIR/SNVs.unfiltered.vcf.gz

# Generate all "interval" files containing those variants.
bedtools intersect -wa -a {intervals_bed} -b $TMPDIR/SVs.mother.vcf.gz \
| sort -k1,1V -k2,2n \
| uniq \
> $TMPDIR/SVs.mother.intervals.bed
bedtools intersect -wa -a {intervals_bed} -b $TMPDIR/SNVs.mother.vcf.gz \
| sort -k1,1V -k2,2n \
| uniq \
> $TMPDIR/SNVs.mother.intervals.bed
bedtools intersect -wa -a {intervals_bed} -b $TMPDIR/SVs.father.vcf.gz \
| sort -k1,1V -k2,2n \
| uniq \
> $TMPDIR/SVs.father.intervals.bed
bedtools intersect -wa -a {intervals_bed} -b $TMPDIR/SNVs.father.vcf.gz \
| sort -k1,1V -k2,2n \
| uniq \
> $TMPDIR/SNVs.father.intervals.bed

# Get the common intervals for both combinations.
comm -12 $TMPDIR/SVs.mother.intervals.bed $TMPDIR/SNVs.father.intervals.bed \
| sort -k1,1V -k2,2n \
| uniq \
> $TMPDIR/SVmo.SNVfa.intervals.bed

comm -12 $TMPDIR/SVs.father.intervals.bed $TMPDIR/SNVs.mother.intervals.bed \
| sort -k1,1V -k2,2n \
| uniq \
> $TMPDIR/SVfa.SNVmo.intervals.bed

# Retrieve the variants corresponding to those intervals.
bedtools intersect -header -wa -a $TMPDIR/SVs.mother.vcf.gz -b $TMPDIR/SVmo.SNVfa.intervals.bed \
| sort -k1,1V -k2,2n \
| uniq \
| bgzip -c \
> $TMPDIR/SVs.mother.intervals.vcf.gz

bedtools intersect -header -wa -a $TMPDIR/SVs.father.vcf.gz -b $TMPDIR/SVfa.SNVmo.intervals.bed \
| sort -k1,1V -k2,2n \
| uniq \
| bgzip -c \
> $TMPDIR/SVs.father.intervals.vcf.gz

bedtools intersect -header -wa -a $TMPDIR/SNVs.mother.vcf.gz -b $TMPDIR/SVfa.SNVmo.intervals.bed \
| sort -k1,1V -k2,2n \
| uniq \
| bgzip -c \
> $TMPDIR/SNVs.mother.intervals.vcf.gz

bedtools intersect -header -wa -a $TMPDIR/SNVs.father.vcf.gz -b $TMPDIR/SVmo.SNVfa.intervals.bed \
| sort -k1,1V -k2,2n \
| uniq \
| bgzip -c \
> $TMPDIR/SNVs.father.intervals.vcf.gz

# Tabix all files
tabix -f $TMPDIR/SVs.mother.intervals.vcf.gz
tabix -f $TMPDIR/SVs.father.intervals.vcf.gz
tabix -f $TMPDIR/SNVs.mother.intervals.vcf.gz
tabix -f $TMPDIR/SNVs.father.intervals.vcf.gz

# Merge (concatenate) them
bcftools concat -a -O z -o {snakemake.output.vcf} \
    $TMPDIR/SVs.mother.intervals.vcf.gz \
    $TMPDIR/SNVs.mother.intervals.vcf.gz \
    $TMPDIR/SVs.father.intervals.vcf.gz \
    $TMPDIR/SNVs.father.intervals.vcf.gz

# Build tabix index.
tabix -f {snakemake.output.vcf}

pushd $(dirname {snakemake.output.vcf})
md5sum $(basename {snakemake.output.vcf}) >$(basename {snakemake.output.vcf}).md5
md5sum $(basename {snakemake.output.tbi}) >$(basename {snakemake.output.tbi}).md5
"""
)
