from math import ceil

from snappy_wrappers.snappy_wrapper import ShellWrapper

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

args = getattr(snakemake.params, "args", {})

# Each of the parallel shards gets an equal part of 75% of the job memory as Java heap
java_mem_mb = ceil(snakemake.resources.mem_mb * 0.75 / args["num_threads"])
annotation_groups = " ".join(f"-G {group}" for group in args["annotation_groups"])

ShellWrapper(snakemake).run(
    r"""
set -x

# Create binning of the reference into windows of roughly the same size.
gatk PreprocessIntervals \
    --reference {snakemake.input.reference} \
    --bin-length {args[window_length]} \
    --output $TMPDIR/raw.interval_list \
    --interval-merging-rule OVERLAPPING_ONLY \
    $(for ignore_contig in {args[ignore_contigs]}; do
        awk "(\$1 ~ /$ignore_contig/) {{ printf(\"--exclude-intervals %s:1-%d\\n\", \$1, \$2) }}" \
            {snakemake.input.reference}.fai
    done)

# Postprocess the Picard-style interval list into properly padded interval strings suitable for
# passing to ``--intervals``.
awk -v PADDING={args[padding]} '
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
run-shard() {{
    job_no=$1
    interval=$2

    GATK_JAVA_MEMORY={java_mem_mb}m
    gatk \
        HaplotypeCaller \
        --java-options "-Xmx$GATK_JAVA_MEMORY -Djava.io.tmpdir=$TMPDIR" \
        --tmp-dir $TMPDIR \
        --output $TMPDIR/shards-output/$(printf %06d $job_no).vcf.gz \
        --reference {snakemake.input.reference} \
        --dbsnp {snakemake.input.dbsnp} \
        {annotation_groups} \
        --intervals $interval \
        {extra_arguments} \
        $(for path in {snakemake.input.bam}; do
            echo --input $path
        done) \
        $(if [[ {args[allow_seq_dict_incompatibility]} == "True" ]]; then
            echo --disable-sequence-dictionary-validation true
        fi)
}}
export -f run-shard

# Perform parallel execution
(
    set -x
    sleep $(echo "scale=3; $RANDOM/32767*10" | bc)s
) # sleep up to 10s to work around bug
num_threads={args[num_threads]}
cat $TMPDIR/final_intervals.txt \
    | parallel -j $num_threads 'run-shard {{#}} {{}}'

# Merge the individual shards' output VCF
bcftools concat \
    --allow-overlaps \
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
