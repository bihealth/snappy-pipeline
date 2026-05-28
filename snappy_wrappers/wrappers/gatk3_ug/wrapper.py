from snappy_wrappers.snappy_wrapper import ShellWrapper

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

args = getattr(snakemake.params, "args", {})

ShellWrapper(snakemake).run(
    r"""
set -x

# Create binning of the reference into windows of roughly the same size.
gatk PreprocessIntervals \
    --reference {snakemake.input.reference} \
    --bin-length {args[window_length]} \
    --output $TMPDIR/raw.interval_list \
    --interval-merging-rule OVERLAPPING_ONLY \
    $(for ignore_chrom in {args[ignore_chroms]}; do \
        awk "(\$1 ~ /$ignore_chrom/) {{ printf(\"--exclude-intervals %s:1-%d\\n\", \$1, \$2) }}" \
            {snakemake.input.reference}.fai; \
    done)

# Postprocess the Picard-style interval list into properly padded interval strings suitable for
# passing to ``--intervals``.
awk -v PADDING=1000 '
    (!/^@/) {{
        chrom=$1;
        start=$2;
        if (start > PADDING) {{
            start = start - PADDING
        }} else {{
            start = 1
        }}
        end=$3;
        printf("%s:%d-%d\n", chrom, start, end);
    }}
' $TMPDIR/raw.interval_list \
> $TMPDIR/final_intervals.txt
wc -l $TMPDIR/final_intervals.txt

# Create per-shard output directory
mkdir -p $TMPDIR/shards-output

# Function to run HaplotypeCaller in "joint" mode - cannot generate gVCF
run-shard()
{{
    job_no=$1
    interval=$2

    GATK_JAVA_MEMORY=3750m
    gatk3 -Xmx$GATK_JAVA_MEMORY -Djava.io.tmpdir=$TMPDIR \
        --analysis_type UnifiedGenotyper \
        --out $TMPDIR/shards-output/$(printf %06d $job_no).vcf.gz \
        --reference_sequence {snakemake.input.reference} \
        --sample_ploidy 2 \
        --dbsnp {snakemake.input.dbsnp} \
        --downsample_to_coverage {args[downsample_to_coverage]} \
        --intervals $interval \
        $(for path in {snakemake.input.bam}; do \
            echo --input_file $path; \
        done) \
        $(if [[ {args[allow_seq_dict_incompatibility]} == "True" ]]; then \
            echo --disable-sequence-dictionary-validation true; \
        fi)
}}
export -f run-shard

# Perform parallel execution
(set -x; sleep $(echo "scale=3; $RANDOM/32767*10" | bc)s) # sleep up to 10s to work around bug
num_threads={args[num_threads]}
cat $TMPDIR/final_intervals.txt \
| parallel --plain -j $num_threads 'run-shard {{#}} {{}}'

# Merge the individual shards' output VCF
bcftools concat \
    --allow-overlaps  \
    -d none \
    -O u \
    $TMPDIR/shards-output/*.vcf.gz \
| bcftools sort \
    -T $TMPDIR/bcftools.sort.XXXXXX \
    -O u \
    /dev/stdin \
| bcftools norm \
    -d exact \
    -f {snakemake.input.reference} \
    -O z \
    -o {snakemake.output.vcf}
tabix {snakemake.output.vcf}
"""
)
