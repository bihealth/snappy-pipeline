"""The ``fgumi/grouping`` wrapper of ``umi_mapping``: from the FASTQ files of a library to its
grouped BAM, fgumi's best-practice chain as one shell pipeline.

``run()`` is the wrapper: it reads the inputs and params, fills the shell texts and runs them
through ``ShellWrapper``.
The job's temp directory (``$TMPDIR``, exported by the wrapper header, ``set -u``) holds the
unmapped BAM, which ``fgumi zipper`` reads twice.
"""

from snappy_wrappers.snappy_wrapper import ShellWrapper
from snappy_wrappers.tools.fgumi import fill, flag, index_prefix

EXTRACT = r"""declare -a reads_left=({reads_left})
declare -a reads_right=({reads_right})
declare -a read_groups=({read_groups})

# Function Definitions ----------------------------------------------------------------------------

# The unmapped BAM of read group $1, with the UMIs moved from the reads into the RX tag
extract_reads() {{
    set -x
    i=$1

    fgumi extract \
        --inputs ${{reads_left[$i]}} ${{reads_right[$i]:-}} \
        --read-structures {read_structures} \
        --read-group-id ${{read_groups[$i]}} \
        --platform-unit ${{read_groups[$i]}} \
        --sample {sample_name} \
        --library {sample_name} \
        --platform {platform} \
        --compression-level 1 \
        --output $TMPDIR/$i.extracted.bam
}}

# The unmapped BAMs of all read groups as one, under one header
merge_read_groups() {{
    set -x

    for f in $TMPDIR/*.extracted.bam; do
        samtools view -H $f > ${{f%.bam}}.hdr.sam
    done
    samtools merge $TMPDIR/merged.hdr.bam $TMPDIR/*.hdr.sam
    samtools cat -h $TMPDIR/merged.hdr.bam -o $TMPDIR/extracted.bam $TMPDIR/*.extracted.bam
    rm $TMPDIR/*.extracted.bam $TMPDIR/*.hdr.sam $TMPDIR/merged.hdr.bam
}}

# Extraction ---------------------------------------------------------------------------------------

for ((i = 0; i < ${{#reads_left[@]}}; i++)); do
    extract_reads $i
done

if [[ ${{#reads_left[@]}} -eq 1 ]]; then
    mv $TMPDIR/0.extracted.bam $TMPDIR/extracted.bam
else
    merge_read_groups
fi"""

CORRECT = r"""# The UMIs corrected against the known UMIs of the kit
fgumi correct \
    --input $TMPDIR/extracted.bam \
    --umi-files {snakemake.input.umis} \
    --max-mismatches {max_mismatches} \
    --min-distance {min_distance} \
    --metrics {snakemake.output.correct_metrics} \
    --compression-level 1 \
    --output $TMPDIR/unmapped.bam
rm $TMPDIR/extracted.bam"""

NO_CORRECT = r"""mv $TMPDIR/extracted.bam $TMPDIR/unmapped.bam"""

GROUP = r"""# Alignment and grouping ---------------------------------------------------------------------------

# fgumi zipper reads the unmapped BAM a second time, so it is a file of the temp directory
fgumi fastq \
    --input $TMPDIR/unmapped.bam \
    --threads {num_threads} \
    | bwa-mem3 mem \
        --bam=0 \
        -p \
        -Y \
        -t {num_threads_align} \
        -K {bases_per_batch} \
        {aligner_extra_arguments} \
        {path_index} \
        - \
    | fgumi zipper \
        --unmapped $TMPDIR/unmapped.bam \
        --reference {snakemake.input.reference} \
        --threads {num_threads} \
        --compression-level 0 \
        --output /dev/stdout \
    | fgumi sort \
        --input /dev/stdin \
        --order template-coordinate \
        --tmp-dir $TMPDIR \
        --threads {num_threads} \
        --compression-level 0 \
        --output - \
    | fgumi group \
        --input /dev/stdin \
        --strategy {strategy} \
        --edits {edits} \
        --min-map-q {min_map_q} \
        --allow-unmapped {allow_unmapped} \
        --realign-window true \
        --metrics {metrics_prefix} \
        --threads {num_threads} \
        --compression-level 1 \
        {group_extra_arguments} \
        --output {snakemake.output.grouped_bam}"""


def command(snakemake) -> str:
    """The shell template of the wrapper: extract (and correct) the UMIs of every read group,
    then align, zipper, sort and group the reads of the library."""
    args = snakemake.params.args
    aligner, group = args["aligner"], args["group"]
    # The reads are inputs. If another task wrote them, the input is that task's out directory and
    # params list the files.
    reads = args.get("input", {})
    input_left = snakemake.input.get("reads_left") or reads["reads_left"]
    input_right = snakemake.input.get("reads_right") or reads.get("reads_right", [])
    correct = args["correct"]
    values = {
        "reads_left": " ".join(input_left),
        "reads_right": " ".join(input_right),
        "read_groups": " ".join(args["read_groups"]),
        "read_structures": " ".join(args["read_structures"]),
        "sample_name": args["sample_name"],
        "platform": args["platform"],
        "num_threads": args["num_threads"],
        "num_threads_align": aligner["num_threads"],
        "bases_per_batch": aligner["bases_per_batch"],
        "aligner_extra_arguments": " ".join(aligner["extra_arguments"]),
        "path_index": index_prefix(snakemake),
        "strategy": args["strategy"],
        "edits": group["edits"],
        "min_map_q": group["min_map_q"],
        "allow_unmapped": flag(group["allow_unmapped"]),
        # fgumi group appends .family_sizes.txt, .grouping_metrics.txt and .position_group_sizes.txt
        "metrics_prefix": snakemake.output.group_family_sizes.removesuffix(".family_sizes.txt"),
        "group_extra_arguments": " ".join(group["extra_arguments"]),
    }
    if correct is not None:
        values["max_mismatches"] = correct["max_mismatches"]
        values["min_distance"] = correct["min_distance"]
    text = "\n\n".join((EXTRACT, CORRECT if correct is not None else NO_CORRECT, GROUP))
    return fill(text, **values)


def run(snakemake) -> None:
    """The wrapper."""
    ShellWrapper(snakemake).run(command(snakemake))
