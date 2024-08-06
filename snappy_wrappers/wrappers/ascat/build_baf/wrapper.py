# -*- coding: utf-8 -*-
"""Wrapper for building BAF files for ASCAT"""

from typing import TYPE_CHECKING

from snakemake.shell import shell

if TYPE_CHECKING:
    from snakemake.script import snakemake

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

shell.executable("/bin/bash")

library_name = getattr(
    snakemake.wildcards,
    "tumor_library_name",
    getattr(snakemake.wildcards, "normal_library_name", None),
)
assert library_name is not None

path_b_af_loci = snakemake.params["args"]["b_af_loci"]
reference_path = snakemake.params["args"]["reference_path"]

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
# Perform pileups at the spot positions.
#
samtools mpileup \
    -l {path_b_af_loci} \
    -I \
    -u \
    -v \
    -t AD \
    -f {reference_path} \
    {snakemake.input.bam} \
| bcftools call \
    -c \
    -O u \
| bcftools query -f "%CHROM\t%POS0\t%END[\t%AD]\n" \
| bgzip -c \
> $TMPDIR/calls.bed.gz
tabix -f $TMPDIR/calls.bed.gz

# -------------------------------------------------------------------------------------------------
# Create VCF file from spots
#
echo "##fileformat=VCFv4.2" \
> $TMPDIR/spots.vcf
echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" \
>> $TMPDIR/spots.vcf

zcat -f {path_b_af_loci} \
| awk \
    -F $'\t' \
    'BEGIN {{ OFS=FS; }}
    ($2 > 0) {{
        print $1, $2 + 1, "SPOT" NR, "N", ".", ".", ".", "."
    }}' \
>> $TMPDIR/spots.vcf
bgzip $TMPDIR/spots.vcf
tabix -f $TMPDIR/spots.vcf.gz

# -------------------------------------------------------------------------------------------------
# Annotate spots VCF and build resulting file from this.
#
cat <<"EOF" >$TMPDIR/ad_header.txt
##INFO=<ID=XAD,Number=1,Type=String,Description="Allelic depth">
EOF

echo -e "\tchrs\tpos\tSAMPLE" \
> {snakemake.output.txt}
bcftools annotate \
    --header-lines $TMPDIR/ad_header.txt \
    --annotations $TMPDIR/calls.bed.gz \
    --columns -CHROM,-FROM,-TO,XAD \
    $TMPDIR/spots.vcf.gz \
| bcftools query \
    -f "%ID\t%CHROM\t%POS\t%XAD\n" \
| sed -e 's/,/\t/g' \
| awk -F $'\t' \
    '
    BEGIN {{ OFS=FS; }}
    (NF == 4) {{ print $1, $2, $3, 0; }}
    (NF == 5 && ($4 + $5 == 0)) {{ print $1, $2, $3, "Inf"; }}
    (NF == 5 && ($4 + $5 > 0)) {{ print $1, $2, $3, $5 / ($4 + $5); }}
    ' \
>> {snakemake.output.txt}

# -------------------------------------------------------------------------------------------------
# Build MD5 sum
#
pushd $(dirname {snakemake.output.txt}) &&
    md5sum $(basename {snakemake.output.txt}) >$(basename {snakemake.output.txt}).md5
"""
)
