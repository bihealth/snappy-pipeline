"""Shell text of the BWA-MEM2 and BWA-MEM3 wrappers, which differ only in the mapper.

A wrapper's template pastes ``read_functions()`` before its mapper function,
which maps the reads of the read group ``$1`` from standard input,
and ``mapping_pipeline()`` after it.
The functions return final shell text: they format the values as Snakemake would.
"""

# The wrappers run with the Python of their environment, so this stays valid Python 3.7.

READ_FUNCTIONS = r"""declare -a reads_left=({reads_left})
declare -a reads_right=({reads_right})
declare -a read_groups=({read_groups})

# Function Definitions ----------------------------------------------------------------------------

# The reads of read group $1, paired reads interleaved with seqtk
#
# In the case of paired reads, write out line counts to temporary files so we can check for problems
# later. A failing zcat in a process substitution does not fail the pipeline, so it leaves a marker.
input_files() {{
    set -x
    i=$1

    if [[ ${{#reads_right[@]}} -eq 0 ]]; then
        zcat ${{reads_left[$i]}}
    else
        seqtk mergepe \
            <({{ zcat ${{reads_left[$i]}} || touch $TMPDIR/zcat_failed; }} | tee >(wc -l > $TMPDIR/count_left.$i)) \
            <({{ zcat ${{reads_right[$i]}} || touch $TMPDIR/zcat_failed; }} | tee >(wc -l > $TMPDIR/count_right.$i))
    fi
}}

# Adapter trimming
trim_adapters() {{
    set -x

    if [[ "{trim_adapters}" == "True" ]]; then
        trimadap-mt -p {num_threads_trimming}
    else
        cat
    fi
}}

# Duplicate masking
mask_duplicates() {{
    set -x

    if [[ "{mask_duplicates}" == "True" ]] && [[ ${{#reads_right[@]}} -gt 0 ]]; then
        samblaster --addMateTags
    elif [[ "{mask_duplicates}" == "True" ]]; then
        # samblaster aborts on reads without mates, unless they are single-end reads by design
        samblaster --addMateTags --ignoreUnmated
    else
        cat
    fi
}}"""

MAPPING_PIPELINE = r"""# The alignments of all read groups, still grouped by read name, under one header
merge_read_groups() {{
    set -x

    for f in $TMPDIR/out.d/*.bam; do
        samtools view -H $f > ${{f%.bam}}.hdr.sam
    done
    samtools merge $TMPDIR/merged.hdr.bam $TMPDIR/out.d/*.hdr.sam
    samtools cat -h $TMPDIR/merged.hdr.bam $TMPDIR/out.d/*.bam \
        | samtools view -h -@ {num_threads_bam_view}
}}

# Alignment postprocessing
#
# This function makes extensive use of tee-piping to prevent the BAM to be read from the disk
# again for statistics.
postproc_bam() {{
    set -x

    samtools sort \
        -T $TMPDIR/sort_bam \
        -m {memory_bam_sort} \
        -@ {num_threads_bam_sort} \
        -O BAM \
        -o /dev/stdout \
        /dev/stdin \
        | tee >(samtools stats /dev/stdin > {output.report_bamstats_txt}) \
        | tee >(samtools flagstat /dev/stdin > {output.report_flagstats_txt}) \
        | tee >(stdbuf -o0 samtools index /dev/stdin -o /dev/stdout > {output.bam_bai}) \
            > {output.bam}

    touch \
        {output.bam} \
        {output.bam_bai}

    # We use stdbuf on samtools index above toether with stdout redirection.  We have seen
    # some interesting timing issues with clustered/distributed file systems and the small
    # index file not having been written by the time that samtools idxstats uses it.  In
    # addition we attempt to run it 5 times with incremental waits to make things more
    # robust.
    for i in {{1..5}}; do
        set +e
        samtools idxstats {output.bam} > {output.report_idxstats_txt}
        ret=$?
        set -e
        if [[ "$ret" -eq 0 ]]; then
            break
        elif [[ "$i" -lt 5 ]]; then
            sleep "${{i}}m"
        fi
    done
    if [[ "$ret" -ne 0 ]]; then
        >&2 echo "samtools idxstats failed 5 times - something wrong with the index?"
        exit 1
    fi
}}

# Check the decompression and the input line counts written out in input_files()
check_input_line_counts() {{
    if [[ -e $TMPDIR/zcat_failed ]]; then
        >&2 echo "Could not decompress the reads."
        exit 1
    fi
    for ((i = 0; i < ${{#reads_left[@]}}; i++)); do
        if [[ -e $TMPDIR/count_left.$i ]]; then
            if ! diff $TMPDIR/count_left.$i $TMPDIR/count_right.$i > /dev/null; then
                >&2 echo "Different number of lines in left/right file of read group ${{read_groups[$i]}}."
                >&2 cat $TMPDIR/count_left.$i
                >&2 cat $TMPDIR/count_right.$i
                exit 1
            fi
        fi
    done
}}

# Run actual tools --------------------------------------------------------------------------------

# Map each read group with its own read group ID; several are concatenated for duplicate masking
if [[ ${{#reads_left[@]}} -eq 1 ]]; then
    input_files 0 \
        | trim_adapters \
        | {mapper} 0 \
        | mask_duplicates \
        | postproc_bam
else
    mkdir -p $TMPDIR/out.d
    for ((i = 0; i < ${{#reads_left[@]}}; i++)); do
        input_files $i \
            | trim_adapters \
            | {mapper} $i \
            | samtools view -b -@ {num_threads_bam_view} -o $TMPDIR/out.d/$(printf %05d $i).bam
    done
    merge_read_groups \
        | mask_duplicates \
        | postproc_bam
fi

check_input_line_counts"""


def _words(value) -> str:
    """``value`` as Snakemake formats it: the items of a sequence separated by spaces."""
    if isinstance(value, (list, tuple, set, frozenset)):
        return " ".join(map(str, value))
    return str(value)


def read_functions(
    reads_left, reads_right, read_groups, trim_adapters, num_threads_trimming, mask_duplicates
) -> str:
    """Declare the reads and read groups as arrays, and define the functions before mapping.

    ``input_files`` writes the reads of a read group, ``trim_adapters`` and ``mask_duplicates``
    filter them.
    """
    return READ_FUNCTIONS.format(
        reads_left=_words(reads_left),
        reads_right=_words(reads_right),
        read_groups=_words(read_groups),
        trim_adapters=trim_adapters,
        num_threads_trimming=num_threads_trimming,
        mask_duplicates=mask_duplicates,
    )


def mapping_pipeline(
    mapper, output, num_threads_bam_view, memory_bam_sort, num_threads_bam_sort
) -> str:
    """Map each read group with the shell function ``mapper``, then sort, index and count.

    ``output`` holds the paths of the BAM, its index and the samtools reports.
    """
    return MAPPING_PIPELINE.format(
        mapper=mapper,
        output=output,
        num_threads_bam_view=num_threads_bam_view,
        memory_bam_sort=memory_bam_sort,
        num_threads_bam_sort=num_threads_bam_sort,
    )
