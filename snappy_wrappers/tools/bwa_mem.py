"""The BWA-MEM2 and BWA-MEM3 wrappers, which differ only in the mapper.

``run()`` is the whole wrapper: it reads the inputs and params, pastes ``read_functions()``
before the mapper function, which maps the reads of the read group ``$1`` from standard input,
and ``mapping_pipeline()`` after it, and runs the script through ``ShellWrapper``.
The text functions return final shell text: they format the values as Snakemake would.
"""

import re

from snappy_wrappers.snappy_wrapper import ShellWrapper

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


#: The template of the wrapper: the read functions, the mapper function and the pipeline
TEMPLATE = """
{read_functions_text}

# Run {mapper_upper} on the reads of read group $1
{mapper_function}() {{
    set -x
    i=$1

    # Decide whether to write split reads as supplementary or secondary (-M means secondary)
    split_as_supp_flag=
    if [[ "{split_as_secondary}" == "True" ]]; then
        split_as_supp_flag="-M"
    fi

    if [[ "{snakemake.params.args[sample_name]}" != "" ]]; then
        rg_arg="-R @RG\\tID:${{read_groups[$i]}}\\tSM:{snakemake.params.args[sample_name]}\\tPL:{snakemake.params.args[platform]}"
    else
        rg_arg=
    fi

    {mapper} mem \\
        {path_index} \\
        $split_as_supp_flag \\
        $rg_arg \\
        -p \\
        -t {num_threads_align} \\
        -K {bases_per_batch} \\
        {extra_arguments} \\
        /dev/stdin
}}

{mapping_pipeline_text}
"""


def _fill(template: str, **values) -> str:
    """``template`` with the placeholders of ``values`` filled in, as Snakemake would format them;
    the other placeholders and the doubled braces stay for ``ShellWrapper``."""
    return re.sub(
        r"\{(" + "|".join(map(re.escape, values)) + r")\}",
        lambda match: str(values[match.group(1)]),
        template,
    )


def run(snakemake, mapper: str) -> None:
    """The wrapper of ``mapper`` (``bwa-mem2`` or ``bwa-mem3``): map the reads of every read
    group of the library, mask duplicates, merge, sort, index and count."""
    args = getattr(snakemake.params, "args", {})

    # The reads are inputs. If another task wrote them, the input is that task's out directory and
    # params list the files.
    reads = args.get("input", {})
    input_left = snakemake.input.get("reads_left") or reads["reads_left"]
    input_right = snakemake.input.get("reads_right") or reads.get("reads_right", "")

    mapper_function = "run_" + mapper.replace("-", "_")
    # The read group IDs, one per read group (file pair) of the library
    read_groups = " ".join(args["read_groups"])
    cmd = _fill(
        TEMPLATE,
        read_functions_text=read_functions(
            input_left,
            input_right,
            read_groups,
            args["trim_adapters"],
            args["num_threads_trimming"],
            args["mask_duplicates"],
        ),
        mapper_upper=mapper.upper(),
        mapper_function=mapper_function,
        mapper=mapper,
        # The index files are inputs, the first one is the .amb file.
        path_index=snakemake.input.index[0].removesuffix(".amb"),
        split_as_secondary=args["split_as_secondary"],
        num_threads_align=args["num_threads_align"],
        bases_per_batch=args["bases_per_batch"],
        mapping_pipeline_text=mapping_pipeline(
            mapper_function,
            snakemake.output,
            args["num_threads_bam_view"],
            args["memory_bam_sort"],
            args["num_threads_bam_sort"],
        ),
    )
    # ShellWrapper fills in the sample name, the platform and the extra arguments
    ShellWrapper(snakemake).run(cmd)
