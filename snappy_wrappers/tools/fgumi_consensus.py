"""The ``fgumi/consensus`` wrapper of ``umi_mapping``: from the grouped BAM of a library to its
filtered, coordinate-sorted consensus BAM, fgumi's best-practice chain as one shell pipeline.

``run()`` is the wrapper: it reads the inputs and params, fills the shell texts and runs them
through ``ShellWrapper``; the tool's ``extra_arguments`` go to the consensus caller.
"""

from snappy_wrappers.snappy_wrapper import ShellWrapper
from snappy_wrappers.tools.fgumi import fill, flag, index_prefix

CLIP = r"""# Function Definitions ----------------------------------------------------------------------------

# The read ends clipped hard, with the per-base consensus tags clipped alike
clip_reads() {{
    set -x

    fgumi clip \
        --input /dev/stdin \
        --reference {snakemake.input.reference} \
        --clipping-mode hard \
        --auto-clip-attributes true \
        --read-one-five-prime {read_one_five_prime} \
        --read-one-three-prime {read_one_three_prime} \
        --read-two-five-prime {read_two_five_prime} \
        --read-two-three-prime {read_two_three_prime} \
        --clip-overlapping-reads {overlapping_reads} \
        --metrics {snakemake.output.clip_metrics} \
        --threads {num_threads} \
        --compression-level 0 \
        {clip_extra_arguments} \
        --output /dev/stdout
}}"""

NO_CLIP = r"""# Function Definitions ----------------------------------------------------------------------------

# No clipping
clip_reads() {{
    cat
}}"""

CONSENSUS = r"""# Consensus calling, realignment and filtering ----------------------------------------------------

# The realigner's bwa-mem3 reads the index next to its reference, so both are linked under one stem
mkdir -p $TMPDIR/ref
ln -s $(realpath {snakemake.input.reference}) $TMPDIR/ref/genome.fa
ln -s $(realpath {snakemake.input.reference_fai}) $TMPDIR/ref/genome.fa.fai
ln -s $(realpath {snakemake.input.reference_dict}) $TMPDIR/ref/genome.dict
for f in {snakemake.input.index}; do
    ln -s $(realpath $f) $TMPDIR/ref/genome.fa${{f#{path_index}}}
done

# realign reverses the per-base tags itself, so filter runs without --reverse-per-base-tags
fgumi {caller} \
    --input {snakemake.input.bam} \
    --min-reads {min_reads} \
    --min-input-base-quality {min_input_base_quality} \
    --output-per-base-tags true \
    --metrics {metrics_prefix} \
    --stats {snakemake.output.consensus_stats} \
    --threads {num_threads} \
    --compression-level 0 \
    {extra_arguments} \
    --output /dev/stdout \
    | fgumi realign \
        --input /dev/stdin \
        --ref $TMPDIR/ref/genome.fa \
        --aligner::command "bwa-mem3 mem --bam=0 -p -Y -t {num_threads_align} -K {chunk_size} {aligner_extra_arguments} {{ref}} /dev/stdin" \
        --aligner::chunk-size {chunk_size} \
        --threads {num_threads} \
        --compression-level 0 \
        --output /dev/stdout \
    | fgumi filter \
        --input /dev/stdin \
        --ref {snakemake.input.reference} \
        --min-reads {filter_min_reads} \
        --max-read-error-rate {max_read_error_rate} \
        --max-base-error-rate {max_base_error_rate} \
        --min-base-quality {min_base_quality} \
        --max-no-call-fraction {max_no_call_fraction} \
        {require_single_strand_agreement} \
        --stats {snakemake.output.filter_stats} \
        --threads {num_threads} \
        --compression-level 0 \
        {filter_extra_arguments} \
        --output /dev/stdout \
    | clip_reads \
    | fgumi sort \
        --input - \
        --order coordinate \
        --write-index true \
        --tmp-dir $TMPDIR \
        --threads {num_threads} \
        --compression-level {compression_level} \
        --output {snakemake.output.bam}

samtools stats {snakemake.output.bam} > {snakemake.output.report_bamstats_txt}
samtools flagstat {snakemake.output.bam} > {snakemake.output.report_flagstats_txt}
samtools idxstats {snakemake.output.bam} > {snakemake.output.report_idxstats_txt}"""


def command(snakemake) -> str:
    """The shell template of the wrapper: call the consensus reads, realign, filter, clip and sort
    them, and run the samtools statistics; ``{extra_arguments}`` go to the caller."""
    args = snakemake.params.args
    aligner, caller, filter_, clip = args["aligner"], args["caller"], args["filter"], args["clip"]
    mode = args["consensus"]
    if mode == "duplex":
        min_reads = ",".join(map(str, caller["min_reads"]))
        filter_min_reads = ",".join(map(str, caller["filter_min_reads"]))
        agreement = "--require-single-strand-agreement {}".format(
            flag(caller["require_single_strand_agreement"])
        )
    else:
        min_reads, filter_min_reads, agreement = caller["min_reads"], caller["filter_min_reads"], ""
    values = {
        "caller": mode,
        "min_reads": min_reads,
        "min_input_base_quality": caller["min_input_base_quality"],
        # The caller appends .family_sizes.txt, .umi_counts.txt and the yield metrics
        "metrics_prefix": snakemake.output.consensus_family_sizes.removesuffix(".family_sizes.txt"),
        "num_threads": args["num_threads"],
        "num_threads_align": aligner["num_threads"],
        "aligner_extra_arguments": " ".join(aligner["extra_arguments"]),
        "path_index": index_prefix(snakemake),
        "chunk_size": args["realign"]["chunk_size"],
        "filter_min_reads": filter_min_reads,
        "max_read_error_rate": filter_["max_read_error_rate"],
        "max_base_error_rate": filter_["max_base_error_rate"],
        "min_base_quality": filter_["min_base_quality"],
        "max_no_call_fraction": filter_["max_no_call_fraction"],
        "require_single_strand_agreement": agreement,
        "filter_extra_arguments": " ".join(filter_["extra_arguments"]),
        "compression_level": args["compression_level"],
    }
    if clip is not None:
        values.update(
            {
                "read_one_five_prime": clip["read_one_five_prime"],
                "read_one_three_prime": clip["read_one_three_prime"],
                "read_two_five_prime": clip["read_two_five_prime"],
                "read_two_three_prime": clip["read_two_three_prime"],
                "overlapping_reads": flag(clip["overlapping_reads"]),
                "clip_extra_arguments": " ".join(clip["extra_arguments"]),
            }
        )
    text = "\n\n".join((CLIP if clip is not None else NO_CLIP, CONSENSUS))
    return fill(text, **values)


def run(snakemake) -> None:
    """The wrapper."""
    ShellWrapper(snakemake).run(command(snakemake))
