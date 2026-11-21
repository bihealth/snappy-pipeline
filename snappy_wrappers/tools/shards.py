"""Shell text of the variant callers that shard their job with GNU parallel.

A wrapper's template pastes ``intervals()`` before its ``run-shard`` function,
``run_shards()`` and ``gather()`` after it.
``run-shard <job number> <interval>`` writes ``$TMPDIR/shards-output/<job number>.vcf.gz``,
with the job number padded to six digits.
The functions return final shell text: they format the values as Snakemake would.
"""

# The wrappers run with the Python of their environment, so this stays valid Python 3.7.

INTERVALS = r"""# Create binning of the reference into windows of roughly the same size.
gatk PreprocessIntervals \
    --reference {reference} \
    --bin-length {window_length} \
    --output $TMPDIR/raw.interval_list \
    --interval-merging-rule OVERLAPPING_ONLY \
    $(for ignore_contig in {ignore_contigs}; do
        awk "(\$1 ~ /$ignore_contig/) {{ printf(\"--exclude-intervals %s:1-%d\\n\", \$1, \$2) }}" \
            {fai}
    done)

# Postprocess the Picard-style interval list into properly padded interval strings suitable for
# passing to ``--intervals``.
awk -v PADDING={padding} '
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
mkdir -p $TMPDIR/shards-output"""

RUN_SHARDS = r"""# Perform parallel execution
(
    set -x
    sleep $(echo "scale=3; $RANDOM/32767*10" | bc)s
) # sleep up to 10s to work around bug
num_threads={num_threads}
cat $TMPDIR/final_intervals.txt \
    | parallel --plain -j $num_threads 'run-shard {{#}} {{}}'"""

GATHER = r"""# Merge the individual shards' output VCF
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
        {check_ref}-f {reference} \
        -O z \
        -o {vcf}
tabix {vcf}"""


def intervals(reference, fai, window_length, padding, ignore_contigs) -> str:
    """Write the intervals of the shards, padded, to ``$TMPDIR/final_intervals.txt``.

    They are windows of ``window_length`` over the reference,
    without the contigs of the ``fai`` that match a regular expression of ``ignore_contigs``.
    """
    return INTERVALS.format(
        reference=reference,
        fai=fai,
        window_length=window_length,
        padding=padding,
        ignore_contigs=" ".join(map(str, ignore_contigs)),
    )


def run_shards(num_threads) -> str:
    """Run ``run-shard`` on each interval, ``num_threads`` at a time."""
    return RUN_SHARDS.format(num_threads=num_threads)


def gather(reference, vcf, check_ref="") -> str:
    """Concatenate, sort and normalize the VCFs of the shards into ``vcf`` and index it.

    ``check_ref`` is the ``-c`` option of ``bcftools norm``.
    """
    # A line of its own before ``-f``
    check_ref = "-c {} \\\n        ".format(check_ref) if check_ref else ""
    return GATHER.format(reference=reference, vcf=vcf, check_ref=check_ref)
