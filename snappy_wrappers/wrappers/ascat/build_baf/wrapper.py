# -*- coding: utf-8 -*-
"""Wrapper for building BAF files for ASCAT"""

from snakemake import shell

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>"

shell.executable("/bin/bash")

library_name = getattr(
    snakemake.wildcards,
    "tumor_library_name",
    getattr(snakemake.wildcards, "normal_library_name", None),
)
assert library_name is not None

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
    -l {snakemake.config[step_config][somatic_purity_ploidy_estimate][ascat][b_af_loci]} \
    -I \
    -u \
    -v \
    -t AD \
    -f {snakemake.config[static_data_config][reference][path]} \
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

zcat -f {snakemake.config[step_config][somatic_purity_ploidy_estimate][ascat][b_af_loci]} \
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
