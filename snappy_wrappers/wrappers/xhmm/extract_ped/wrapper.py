# -*- coding: utf-8 -*-
# Extract pedigree members only from VCF file.

from snakemake.shell import shell

shell(
    r"""
# -----------------------------------------------------------------------------
# Redirect stderr to log file by default and enable printing executed commands
exec &> >(tee -a "{snakemake.log}")
set -x
# -----------------------------------------------------------------------------

set -euo pipefail

export TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT

echo '{snakemake.params.ped_members}' \
| tr ' ' '\n' \
| LANG=C sort \
>$TMPDIR/samples_raw.txt

sort {snakemake.input.filtered_samples} \
| LANG=C sort \
>$TMPDIR/filtered_samples.txt

cat >$TMPDIR/header.txt <<"EOF"
##INFO=<ID=SVMETHOD,Number=1,Type=String,Description="Type of approach used to detect SV">
EOF

# We need to exclude the filtered samples from extracting PED.  For the record, we will ad
# a header adding a note bout this sample, though.
comm -23 $TMPDIR/samples_raw.txt $TMPDIR/filtered_samples.txt \
> $TMPDIR/samples.txt

if ! diff $TMPDIR/samples_raw.txt $TMPDIR/samples.txt; then
    comm -23 $TMPDIR/samples_raw.txt $TMPDIR/samples.txt \
    | awk '{{ print "##SAMPLE=<ID=" $1 ",Warning=\"Filtered out by XHMM\">"}}' \
    >>$TMPDIR/header.txt
fi

head -n 1000 $TMPDIR/samples_raw.txt $TMPDIR/filtered_samples.txt $TMPDIR/samples.txt $TMPDIR/header.txt

if grep . $TMPDIR/samples.txt; then
    bcftools view \
        --force-samples \
        --samples-file $TMPDIR/samples.txt \
        --output-type u \
        {snakemake.input.vcf} \
    | bcftools annotate \
        --header-lines $TMPDIR/header.txt \
    | awk -F $'\t' '
        BEGIN {{ OFS = FS; }}
        /^#/ {{ print $0; }}
        /^[^#]/ {{ $8 = $8 ";SVMETHOD=XHMMv0.0.0.2016_01_04.cc14e52"; print $0; }} ' \
    | bcftools view \
        --output-file {snakemake.output.vcf} \
        --output-type z \
        --include '(GT == "alt")'
else
    # no sample, just keep header
    bcftools view \
        --force-samples \
        --samples-file $TMPDIR/samples.txt \
        {snakemake.input.vcf} \
    | zgrep '^#' \
    | bgzip -c \
    > {snakemake.output.vcf}
fi

tabix -f {snakemake.output.vcf}

pushd $(dirname {snakemake.output.vcf})
md5sum $(basename {snakemake.output.vcf}) >$(basename {snakemake.output.vcf}).md5
md5sum $(basename {snakemake.output.vcf}).tbi >$(basename {snakemake.output.vcf}).tbi.md5
"""
)
