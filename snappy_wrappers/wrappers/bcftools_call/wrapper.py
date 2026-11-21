from typing import TYPE_CHECKING

from snappy_wrappers.snappy_wrapper import ShellWrapper
from snappy_wrappers.tools import shards

if TYPE_CHECKING:
    from snakemake.iocontainers import snakemake

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"
args = getattr(snakemake.params, "args", {})

reference_path = snakemake.input.reference
reference_index_path = snakemake.input.reference_fai

assembly = args["assembly"]

ignore_contigs = args["ignore_contigs"]

max_depth = args["max_depth"]
max_indel_depth = args["max_indel_depth"]

window_length = args["window_length"]
padding = args["padding"]
annotate = args["annotate"]
per_sample_mf = "--per-sample-mF" if args["per_sample_mf"] else ""
redo_baq = "--redo-BAQ" if args["redo_baq"] else ""
variants_only = "--variants-only" if args["variants_only"] else ""
num_threads = args["num_threads"]

shard_intervals = shards.intervals(
    reference_path, reference_index_path, window_length, padding, ignore_contigs
)
run_shards = shards.run_shards(num_threads)
gather_shards = shards.gather(reference_path, snakemake.output.vcf)

ShellWrapper(snakemake).run(
    r"""
# Run actual tools --------------------------------------------------------------------------------

{shard_intervals}

# Function to run "bcftools call"
run-shard() {{
    set -x

    job_no=$1
    interval=$2

    bcftools mpileup \
        -Ou \
        --annotate {annotate} \
        -f {reference_path} \
        {per_sample_mf} \
        --max-depth {max_depth} \
        --max-idepth {max_indel_depth} \
        {redo_baq} \
        --regions $2 \
        {snakemake.input.bam} \
        | bcftools call \
            --multiallelic-caller \
            $(if [[ "{assembly}" != "unknown" ]]; then
                echo --ploidy={assembly}
                echo --samples-file {snakemake.input.ped}
            fi) \
            {variants_only} \
            {extra_arguments} \
            -Oz \
            -o $TMPDIR/shards-output/$(printf %06d $job_no).vcf.gz
    tabix -f $TMPDIR/shards-output/$(printf %06d $job_no).vcf.gz
}}
export -f run-shard

{run_shards}

{gather_shards}
"""
)
