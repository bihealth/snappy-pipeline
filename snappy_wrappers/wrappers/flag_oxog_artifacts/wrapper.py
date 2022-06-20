# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for additional flagging.

TODO: rename to reflect that mroe is done than just flagging the oxog artifacts.
"""

from snakemake import shell

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

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

# Extract genotype calls of the somatic DNA library and pipe into the
# actual filtration command. We simply build a BED file which we then use
# below for flagging.
bcftools view \
    -s "{snakemake.wildcards.tumor_library}" \
    -O u \
    -i '
        (AF[0] <= 0.08 || AD[1] < 5) &&
            ((REF == "G" && ALT ~ "T") || (REF = "C" && ALT ~ "A"))' \
    {snakemake.input.vcf} \
| bcftools query \
    -i 'GT ~ "1"' \
    -f "%CHROM\t%POS0\t%END\t1\n" \
| bgzip -c \
> $TMPDIR/oxog_positions.bed.gz
tabix -f $TMPDIR/oxog_positions.bed.gz

# Create BED file with maximal alternate allele depth of normal sample.
bcftools query -f "%CHROM\t%POS0\t%END\t[%AD\t]\n" {snakemake.input.vcf} \
| awk -F $'\t' '
    BEGIN {{
        OFS = FS;
    }}
    {{
        split($5, arr, ",");
        n = 0;
        for (i = 2; i <= length(arr); ++i) {{
            if (arr[i] > n) {{
                n = arr[i];
            }}
    }}
    print $1, $2, $3, n;
}}' \
| bgzip -c \
> $TMPDIR/alt_mdaa.bed.gz
tabix -f $TMPDIR/alt_mdaa.bed.gz

# Build header lines to add
cat <<"EOF" >$TMPDIR/header.txt
##FORMAT=<ID=CUBI_OXOG,Number=1,Type=Integer,Description="Simple heuristic for oxo-G filter (low support thresholds)">
##INFO=<ID=CUBI_NMDAA,Number=1,Type=Integer,Description="Maximal depth of all alternate alleles in normal sample">
EOF

# Flag the input VCF file with the BED files we built earlier.
bcftools annotate \
    -s "{snakemake.wildcards.tumor_library}" \
    --header-lines $TMPDIR/header.txt \
    --annotations $TMPDIR/oxog_positions.bed.gz \
    --columns "-CHROM,-FROM,-TO,FMT/CUBI_OXOG" \
    -O u \
    {snakemake.input.vcf} \
| bcftools annotate \
    --annotations $TMPDIR/alt_mdaa.bed.gz \
    --columns "-CHROM,-FROM,-TO,CUBI_NMDAA" \
    -O z \
    -o {snakemake.output.vcf}
tabix -f {snakemake.output.vcf}

pushd $(dirname {snakemake.output.vcf}) && \
    md5sum $(basename {snakemake.output.vcf}) >$(basename {snakemake.output.vcf}).md5 && \
    md5sum $(basename {snakemake.output.tbi}) >$(basename {snakemake.output.tbi}).md5 && \
    popd
"""
)
