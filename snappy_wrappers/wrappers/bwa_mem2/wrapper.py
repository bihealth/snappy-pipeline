from typing import TYPE_CHECKING

from snappy_wrappers.snappy_wrapper import ShellWrapper

if TYPE_CHECKING:
    from snakemake.iocontainers import snakemake

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

args = getattr(snakemake.params, "args", {})

# The reads are inputs. If another task wrote them, the input is that task's out directory and
# params list the files.
reads = args.get("input", {})
input_left = snakemake.input.get("reads_left") or reads["reads_left"]
input_right = snakemake.input.get("reads_right") or reads.get("reads_right", "")

# The index files are inputs, the first one is the .amb file.
path_index = snakemake.input.index[0].removesuffix(".amb")
trim_adapters = args["trim_adapters"]
num_threads_trimming = args["num_threads_trimming"]
mask_duplicates = args["mask_duplicates"]
num_threads_bam_view = args["num_threads_bam_view"]
memory_bam_sort = args["memory_bam_sort"]
num_threads_bam_sort = args["num_threads_bam_sort"]
num_threads_align = args["num_threads_align"]
bases_per_batch = args["bases_per_batch"]
split_as_secondary = args["split_as_secondary"]
# The read group IDs, one per read group (file pair) of the library
read_groups = " ".join(args["read_groups"])

ShellWrapper(snakemake).run(
    r"""
set -x

declare -a reads_left=({input_left})
declare -a reads_right=({input_right})
declare -a read_groups=({read_groups})

# Function Definitions ----------------------------------------------------------------------------

# The reads of read group $1, paired reads interleaved with seqtk
#
# In the case of paired reads, write out line counts to temporary files so we can check for problems
# later. A failing zcat in a process substitution does not fail the pipeline, so it leaves a marker.
input_files()
{{
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
trim_adapters()
{{
    set -x

    if [[ "{trim_adapters}" == "True" ]]; then
        trimadap-mt -p {num_threads_trimming}
    else
        cat
    fi
}}

# Duplicate masking
mask_duplicates()
{{
    set -x

    if [[ "{mask_duplicates}" == "True" ]]; then
        samblaster --addMateTags
    else
        cat
    fi
}}

# Run BWA-MEM2 on the reads of read group $1
run_bwa_mem2()
{{
    set -x
    i=$1

    # Decide whether to write split reads as supplementary or secondary (-M means secondary)
    split_as_supp_flag=
    if [[ "{split_as_secondary}" == "True" ]]; then
        split_as_supp_flag="-M"
    fi

    if [[ "{snakemake.params.args[sample_name]}" != "" ]]; then
        rg_arg="-R @RG\tID:${{read_groups[$i]}}\tSM:{snakemake.params.args[sample_name]}\tPL:{snakemake.params.args[platform]}"
    else
        rg_arg=
    fi

    bwa-mem2 mem \
        {path_index} \
        $split_as_supp_flag \
        $rg_arg \
        -p \
        -t {num_threads_align} \
        -K {bases_per_batch} \
        {extra_arguments} \
        /dev/stdin
}}

# The alignments of all read groups, still grouped by read name, under one header
merge_read_groups()
{{
    set -x

    for f in $TMPDIR/out.d/*.bam; do
        samtools view -H $f >${{f%.bam}}.hdr.sam
    done
    samtools merge $TMPDIR/merged.hdr.bam $TMPDIR/out.d/*.hdr.sam
    samtools cat -h $TMPDIR/merged.hdr.bam $TMPDIR/out.d/*.bam \
    | samtools view -h -@ {num_threads_bam_view}
}}

# Alignment postprocessing
#
# This function makes extensive use of tee-piping to prevent the BAM to be read from the disk
# again for statistics.
postproc_bam()
{{
    set -x

    samtools sort \
        -T $TMPDIR/sort_bam \
        -m {memory_bam_sort} \
        -@ {num_threads_bam_sort} \
        -O BAM \
        -o /dev/stdout \
        /dev/stdin \
    | tee >(samtools stats    /dev/stdin >{snakemake.output.report_bamstats_txt}) \
    | tee >(samtools flagstat /dev/stdin >{snakemake.output.report_flagstats_txt}) \
    | tee >(stdbuf -o0 samtools index /dev/stdin -o /dev/stdout >{snakemake.output.bam_bai}) \
    > {snakemake.output.bam}

    touch \
        {snakemake.output.bam} \
        {snakemake.output.bam_bai}

    # We use stdbuf on samtools index above toether with stdout redirection.  We have seen
    # some interesting timing issues with clustered/distributed file systems and the small
    # index file not having been written by the time that samtools idxstats uses it.  In
    # addition we attempt to run it 5 times with incremental waits to make things more
    # robust.
    for i in {{1..5}}; do
        set +e
        samtools idxstats {snakemake.output.bam} >{snakemake.output.report_idxstats_txt}
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
check_input_line_counts()
{{
    if [[ -e $TMPDIR/zcat_failed ]]; then
        >&2 echo "Could not decompress the reads."
        exit 1
    fi
    for ((i = 0; i < ${{#reads_left[@]}}; i++)); do
        if [[ -e $TMPDIR/count_left.$i ]]; then
            if ! diff $TMPDIR/count_left.$i $TMPDIR/count_right.$i >/dev/null; then
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
    | run_bwa_mem2 0 \
    | mask_duplicates \
    | postproc_bam
else
    mkdir -p $TMPDIR/out.d
    for ((i = 0; i < ${{#reads_left[@]}}; i++)); do
        input_files $i \
        | trim_adapters \
        | run_bwa_mem2 $i \
        | samtools view -b -@ {num_threads_bam_view} -o $TMPDIR/out.d/$(printf %05d $i).bam
    done
    merge_read_groups \
    | mask_duplicates \
    | postproc_bam
fi

check_input_line_counts
"""
)
