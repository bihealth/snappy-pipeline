# -*- coding: utf-8 -*-
"""Wrapper for building CNV files for ASCAT"""

from snakemake import shell

__author__ = "Clemens Messerschmidt <clemens.messerschmidt@bih-charite.de>"

shell.executable("/bin/bash")

library_name = getattr(
    snakemake.wildcards,
    "tumor_library_name",
    getattr(snakemake.wildcards, "normal_library_name", None),
)
assert library_name is not None

# tumor log2 read counts are in columns 5 and 6 in the log2_read_counts.igv file
if snakemake.wildcards.get("tumor_library_name", False):
    log2_column = 5
else:
    log2_column = 6

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

zcat -f {snakemake.config[step_config][somatic_purity_ploidy_estimate][ascat][b_af_loci]} \
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
# Build BED file for annotating SNP positions
#
bins={snakemake.input.bins}
log=$(dirname $bins)/../CNAprofiles/log2_read_counts.igv

tail -n +3 ${{log}} \
| awk -F $'\t' 'BEGIN {{OFS = FS}} {{print $1, $2-1, $3, ${log2_column} }}' \
| bgzip -c \
> $TMPDIR/cov0.bed.gz
tabix -f $TMPDIR/cov0.bed.gz

# -------------------------------------------------------------------------------------------------
# Annotate spots VCF and build resulting file from this.
#
cat <<"EOF" >$TMPDIR/ad_header.txt
##INFO=<ID=NCOV0,Number=1,Type=String,Description="Normalized coverage (incl. q0 reads)">
EOF

echo -e "\tchrs\tpos\tSAMPLE" \
> {snakemake.output.txt}
bcftools annotate \
    --header-lines $TMPDIR/ad_header.txt \
    --annotations $TMPDIR/cov0.bed.gz \
    --columns -CHROM,-FROM,-TO,NCOV0 \
    $TMPDIR/spots.vcf.gz \
| bcftools query \
    -f "%ID\t%CHROM\t%POS\t%NCOV0\n" \
| sed -e 's/,/\t/g' \
| awk -F $'\t' '
    BEGIN {{ OFS=FS; }}
    ($2 == 0 || $4 == ".") {{ next; }}
    {{ print $1, $2, $3, $4; }}
    ' \
>> {snakemake.output.txt}

# -------------------------------------------------------------------------------------------------
# Build MD5 sum
#
pushd $(dirname {snakemake.output.txt}) &&
    md5sum $(basename {snakemake.output.txt}) >$(basename {snakemake.output.txt}).md5
"""
)
