# -*- coding: utf-8 -*-
"""Wrapper for running GATK ReadBackedPhasing in parallel, genome is split into windows
"""

from snakemake.shell import shell

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

shell.executable("/bin/bash")

shell(
    r"""
set -x

# Hack: get back bin directory of base/root environment.
export PATH=$PATH:$(dirname $(dirname $(which conda)))/bin

# Allow this number of variant per 10k window.
PER10K=100

# Also pipe everything to log file
if [[ -n "{snakemake.log}" ]]; then
    if [[ "$(set +e; tty; set -e)" != "" ]]; then
        rm -f "{snakemake.log}" && mkdir -p $(dirname {snakemake.log})
        exec &> >(tee -a "{snakemake.log}" >&2)
    else
        rm -f "{snakemake.log}" && mkdir -p $(dirname {snakemake.log})
        echo "No tty, logging disabled" >"{snakemake.log}"
    fi
fi

# Create auto-cleaned temporary directory
export TMPDIR=$(mktemp -d)
#trap "rm -rf $TMPDIR" EXIT

# Build list of intervals with <1% variants -------------------------------------------------------

# Generate variant counts of 10kbp intervals
bcftools query \
    -r $(echo {snakemake.params[args][intervals]} | tr -d ',' | tr ' ' ',') \
    -f "%CHROM\t%POS\n" \
    {snakemake.input.vcf} \
| awk -F $'\t' '
        BEGIN {{ OFS = FS; prev_chrom = 0; prev_count = 0; }}
        {{
            $2 -= 1;
            if ($1 == prev_chrom &&
                int($2 / 10000) == prev_window) {{
                prev_count += 1;
            }} else {{
                if (prev_chrom != 0) {{
                    print prev_chrom, (prev_window * 10000),
                            (prev_window + 1) * 10000, prev_count;
                }}
                prev_chrom = $1;
                prev_window = int($2 / 10000);
                prev_count = 0;
            }}
        }}
        END {{
                if (prev_chrom != 0) {{
                    print prev_chrom, (prev_window * 10000),
                            (prev_window + 1) * 10000, prev_count;
                }}
        }}' \
> $TMPDIR/var_counts.bed

# Cut out intervals to perform phasing for
awk \
    -v per10k=$PER10K \
    -F $'\t' \
    'BEGIN {{ OFS = FS; }} ($4 > per10k) {{ print }}' \
    $TMPDIR/var_counts.bed \
| bedtools merge \
> $TMPDIR/intervals.neg.bed

# Cut out intervals to NOT perform phasing for
awk \
    -v per10k=$PER10K \
    -F $'\t' \
    'BEGIN {{ OFS = FS; }} ($4 <= per10k) {{ print }}' \
    $TMPDIR/var_counts.bed \
| bedtools merge \
> $TMPDIR/intervals.pos.bed

head -n 1000 $TMPDIR/intervals.*.bed

ls -lh $TMPDIR

if [[ -s $TMPDIR/intervals.neg.bed ]]; then
    arg="-R $TMPDIR/intervals.neg.bed"
else
    arg="-r ''"
fi

bcftools view $arg {snakemake.input.vcf} \
| bgzip -c > $TMPDIR/intervals.neg.vcf.gz
tabix $TMPDIR/intervals.neg.vcf.gz

bcftools view -H $TMPDIR/intervals.neg.vcf.gz | wc -l

if [[  -s $TMPDIR/intervals.pos.bed ]]; then
    arg="-R $TMPDIR/intervals.pos.bed"
else
    arg="-r ''"
fi

bcftools view $arg {snakemake.input.vcf} \
| snappy-vcf_filter_to_info \
    --input-vcf /dev/stdin \
    --output-vcf /dev/stdout \
| bgzip -c > $TMPDIR/intervals.pos.nofilter.vcf.gz
tabix $TMPDIR/intervals.pos.nofilter.vcf.gz

bcftools view -H $TMPDIR/intervals.pos.nofilter.vcf.gz | wc -l

# Call Wrapper ------------------------------------------------------------------------------------

export MALLOC_ARENA_MAX=4

gatk_nonfree \
    -Xmx8g \
    -Djava.io.tmpdir=$TMPDIR \
    --analysis_type ReadBackedPhasing \
    -nct 1 \
    --phaseQualityThresh {snakemake.config[step_config][variant_phasing][gatk_read_backed_phasing][phase_quality_threshold]} \
    --variant $TMPDIR/intervals.pos.nofilter.vcf.gz \
    --out $TMPDIR/phased.nofilter.vcf.gz \
    $(for bam in {snakemake.input.bam}; do echo -I $bam; done) \
    --reference_sequence {snakemake.config[static_data_config][reference][path]} \
    -L $(echo {snakemake.params[args][intervals]} | tr -d ',' | tr ' ' ',')

snappy-vcf_filter_from_info \
    --input-vcf $TMPDIR/phased.nofilter.vcf.gz \
    --output-vcf $TMPDIR/phased.vcf.gz
tabix -f $TMPDIR/phased.vcf.gz

# Merge Resulting Variants ------------------------------------------------------------------------

bcftools \
    concat \
    -D \
    -a \
    $TMPDIR/phased.vcf.gz \
    $TMPDIR/intervals.neg.vcf.gz \
    -O z \
    -o {snakemake.output.vcf}

tabix -f {snakemake.output.vcf}

pushd $(dirname {snakemake.output.vcf}) && \
    md5sum $(basename {snakemake.output.vcf}) >$(basename {snakemake.output.vcf}).md5 && \
    md5sum $(basename {snakemake.output.vcf_tbi}) >$(basename {snakemake.output.vcf_tbi}).md5
"""
)
