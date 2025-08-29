# -*- coding: utf-8 -*-
"""Wrapper for building CNV files for ASCAT"""

from typing import TYPE_CHECKING

from snakemake.shell import shell

if TYPE_CHECKING:
    from snakemake.script import snakemake

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

shell.executable("/bin/bash")

args = getattr(snakemake.params, "args", {})

shell(
    r"""
set -x

export TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT

# Also pipe stderr to log file
if [[ -n "{snakemake.log}" ]]; then
    if [[ "$(set +e; tty; set -e)" != "" ]]; then
        rm -f "{snakemake.log}" && mkdir -p $(dirname {snakemake.log})
        exec 2> >(tee -a "{snakemake.log}" >&2)
    else
        rm -f "{snakemake.log}" && mkdir -p $(dirname {snakemake.log})
        echo "No tty, logging disabled" >"{snakemake.log}"
    fi
fi

# -------------------------------------------------------------------------------------------------
# Create VCF file from spots
#
echo "##fileformat=VCFv4.2" \
> $TMPDIR/spots.vcf
echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" \
>> $TMPDIR/spots.vcf

zcat -f {args[path_b_af_loci]} \
| awk \
    -F $'\t' '
    BEGIN {{ OFS=FS; }}
    ($2 > 0) {{
        print $1, $2 + 1, "SPOT" NR, "N", ".", ".", ".", "."
    }}' \
>> $TMPDIR/spots.vcf
bgzip $TMPDIR/spots.vcf
tabix -f $TMPDIR/spots.vcf.gz

# -------------------------------------------------------------------------------------------------
# Build temporary BCF file with the normalized coverage information.
#
# TODO: should become a conda package!
/fast/groups/cubi/scratch/mholtgr/cnvetti quick wgs-cov-bins \
    --reference {args[reference]} \
    --input {snakemake.input.bam} \
    --output $TMPDIR/cov.bcf

# -------------------------------------------------------------------------------------------------
# Build BED file for annotating, does not work with BCF :(
#
bcftools query \
    -f "%CHROM\t%POS0\t%END[\t%NCOV]\n" \
    $TMPDIR/cov.bcf \
| bgzip -c \
> $TMPDIR/ncov.bed.gz
tabix -f $TMPDIR/ncov.bed.gz

# -------------------------------------------------------------------------------------------------
# Annotate spots VCF and build resulting file from this.
#
cat <<"EOF" >$TMPDIR/ad_header.txt
##INFO=<ID=NCOV,Number=1,Type=String,Description="Normalized coverage (incl. q0 reads)">
EOF

echo -e "\tchrs\tpos\tSAMPLE" \
> {snakemake.output.txt}
bcftools annotate \
    --header-lines $TMPDIR/ad_header.txt \
    --annotations $TMPDIR/ncov.bed.gz \
    --columns -CHROM,-FROM,-TO,NCOV \
    $TMPDIR/spots.vcf.gz \
| bcftools query \
    -f "%ID\t%CHROM\t%POS\t%NCOV\n" \
| sed -e 's/,/\t/g' \
| awk -F $'\t' '
    BEGIN {{ OFS=FS; }}
    ($2 == 0 ) {{ next; }}
    {{ print $1, $2, $3, log($4) / log(2); }}
    ' \
>> {snakemake.output.txt}

# -------------------------------------------------------------------------------------------------
# Build MD5 sum
#
pushd $(dirname {snakemake.output.txt}) &&
    md5sum $(basename {snakemake.output.txt}) >$(basename {snakemake.output.txt}).md5
"""
)
