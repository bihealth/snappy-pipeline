# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for het_comp filter for variant_filtration."""

import os
import sys

from snakemake import shell

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

# Prelude -----------------------------------------------------------------------------------------

shell.executable("/bin/bash")
shell.prefix("set -eu -o pipefail -x; ")

# Get path to this file's (wrapper.py) directory.
base_dir = os.path.dirname(os.path.realpath(__file__))

args = getattr(snakemake.params, "args", {})

# Short-circuit in case of performing no filtration
if args["filter_mode"] == "passthrough":
    shell(
        r"""
    # Het. comp. mode set to "passthrough", just copy out the data.
    ln -sr {snakemake.input.vcf} {snakemake.output.vcf}
    ln -sr {snakemake.input.vcf_md5} {snakemake.output.vcf_md5}
    ln -sr {snakemake.input.vcf_tbi} {snakemake.output.vcf_tbi}
    ln -sr {snakemake.input.vcf_tbi_md5} {snakemake.output.vcf_tbi_md5}
    """
    )
    sys.exit(0)  # everything went well!


# Actual Filtration -------------------------------------------------------------------------------

shell(
    r"""
set -x

# Setup auto-cleaned TMPDIR
export TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT

# Load library with helper functions.
source {base_dir}/../../wgs_sv_filtration/funcs.sh

# Get name and number of index, father, and mother.
index={args[index_library]}
father=$(awk '($2 == "'$index'") {{ print $3; }}' {snakemake.input.ped})
mother=$(awk '($2 == "'$index'") {{ print $4; }}' {snakemake.input.ped})

index_no=$(get_index {snakemake.input.vcf} "$index")
father_no=$(get_index {snakemake.input.vcf} "$father")
mother_no=$(get_index {snakemake.input.vcf} "$mother")

## XXX perform actual filtration XXX

### Missing parents should not exist or? Giving the inheritance filter I used.

### Create par1 file

# Build filter string for inclusion.
include="(GT[$index_no] == \"0/1\")"
test -n "$father_no" && include+=" && (GT[$father_no] == \"0/1\")"
test -n "$mother_no" && include+=" && (GT[$mother_no] == \"0/0\")"
# Perform the actual filtration.
bcftools view \
    -i "$include" \
    -O z \
    -o $TMPDIR/par1.SNV.vcf.gz \
    {snakemake.input.vcf}

### Create par2 file

# Build filter for inclusion.
include="(GT[$index_no] == \"0/1\")"
test -n "$father_no" && include+=" && (GT[$father_no] == \"0/0\")"
test -n "$mother_no" && include+=" && (GT[$mother_no] == \"0/1\")"
# Perform the actual filtration.
bcftools view \
    -i "$include" \
    -O z \
    -o $TMPDIR/par2.SNV.vcf.gz \
    {snakemake.input.vcf}

### Determine intervals to use for hetcomp criteria

if [[ "{args[filter_mode]}" == tads ]]; then
    intervals_bed={args[filter_config][all_tads]}

elif [[ "{args[filter_mode]}" == intervals500 ]]; then
    ### is it ok to use any parent?
    zcat $TMPDIR/par1.SNV.vcf.gz \
    | {{ grep -v ^\# || true; }} \
    | awk -F'\t' 'BEGIN {{ OFS = FS }} {{ left = 500; if (left > $2) {{ left = $2 }} print $1, $2 - left, $2 + 500 }}' \
    > $TMPDIR/par1.SNV.bed
    intervals_bed=$TMPDIR/par1.SNV.bed
elif [[ "{args[filter_mode]}" == gene ]]; then
    intervals_bed={args[filter_config][all_genes]}
fi

echo $intervals_bed

### Add exonic + effect filter for ARHC in genes

if [[ "{args[filter_mode]}" == gene ]]; then
    bcftools filter \
        -i '(INFO/ANN ~ "MODERATE") || (INFO/ANN ~ "HIGH")' \
        -O z \
        -o $TMPDIR/par1.SNV.toUse.vcf.gz \
        $TMPDIR/par1.SNV.vcf.gz
    bcftools filter \
        -i '(INFO/ANN ~ "MODERATE") || (INFO/ANN ~ "HIGH")' \
        -O z \
        -o $TMPDIR/par2.SNV.toUse.vcf.gz \
        $TMPDIR/par2.SNV.vcf.gz

else
    mv $TMPDIR/par1.SNV.vcf.gz $TMPDIR/par1.SNV.toUse.vcf.gz
    mv $TMPDIR/par2.SNV.vcf.gz $TMPDIR/par2.SNV.toUse.vcf.gz
fi

## XXX perform actual filtration XXX

# Generate list of affected regions per parent
# NB: this stupid 8fir bedfiles) sort command has to be used because comm requires it

bedtools intersect -wa -a $intervals_bed -b $TMPDIR/par1.SNV.toUse.vcf.gz \
| sort -k1,1V -k2,2n \
| uniq \
> $TMPDIR/par1.OKregions.bed
bedtools intersect -wa -a $intervals_bed -b $TMPDIR/par2.SNV.toUse.vcf.gz \
| sort -k1,1V -k2,2n \
| uniq \
> $TMPDIR/par2.OKregions.bed

# Find common regions
# NB: we then resort on natural sorting for bed files

comm -12 $TMPDIR/par1.OKregions.bed $TMPDIR/par2.OKregions.bed \
| sort -k1,1V -k2,2n \
| uniq \
> $TMPDIR/common.OKregions.bed

# Filter the variant files
bedtools intersect -header -wa -a $TMPDIR/par1.SNV.toUse.vcf.gz -b $TMPDIR/common.OKregions.bed \
| sort -k1,1V -k2,2n \
| uniq \
| bgzip -c > $TMPDIR/par1.SNV.ARHC.vcf.gz
bedtools intersect -header -wa -a $TMPDIR/par2.SNV.toUse.vcf.gz -b $TMPDIR/common.OKregions.bed \
| sort -k1,1V -k2,2n \
| uniq \
| bgzip -c > $TMPDIR/par2.SNV.ARHC.vcf.gz

tabix -f $TMPDIR/par1.SNV.ARHC.vcf.gz
tabix -f $TMPDIR/par2.SNV.ARHC.vcf.gz

# Concatenate the variant files
bcftools concat -a -O z -o {snakemake.output.vcf} \
$TMPDIR/par1.SNV.ARHC.vcf.gz $TMPDIR/par2.SNV.ARHC.vcf.gz

tabix -f {snakemake.output.vcf}

pushd $(dirname {snakemake.output.vcf})
md5sum $(basename {snakemake.output.vcf}) >$(basename {snakemake.output.vcf}).md5
md5sum $(basename {snakemake.output.vcf_tbi}) >$(basename {snakemake.output.vcf_tbi}).md5
"""
)
