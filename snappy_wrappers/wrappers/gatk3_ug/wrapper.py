from snappy_wrappers.snappy_wrapper import ShellWrapper
from snappy_wrappers.tools import shards
from snappy_wrappers.tools.java import java_heap_mb

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

args = getattr(snakemake.params, "args", {})

# Each of the parallel shards gets an equal part of 75% of the job memory as Java heap
java_mem_mb = java_heap_mb(snakemake, shards=args["num_threads"])
shard_intervals = shards.intervals(
    snakemake.input.reference,
    f"{snakemake.input.reference}.fai",
    args["window_length"],
    args["padding"],
    args["ignore_contigs"],
)
run_shards = shards.run_shards(args["num_threads"])
gather_shards = shards.gather(snakemake.input.reference, snakemake.output.vcf)

ShellWrapper(snakemake).run(
    r"""
{shard_intervals}

# Function to run HaplotypeCaller in "joint" mode - cannot generate gVCF
run-shard() {{
    job_no=$1
    interval=$2

    GATK_JAVA_MEMORY={java_mem_mb}m
    gatk3 -Xmx$GATK_JAVA_MEMORY -Djava.io.tmpdir=$TMPDIR \
        --analysis_type UnifiedGenotyper \
        --out $TMPDIR/shards-output/$(printf %06d $job_no).vcf.gz \
        --reference_sequence {snakemake.input.reference} \
        --sample_ploidy {args[sample_ploidy]} \
        --dbsnp {snakemake.input.dbsnp} \
        --downsample_to_coverage {args[downsample_to_coverage]} \
        --intervals $interval \
        {extra_arguments} \
        $(for path in {snakemake.input.bam}; do
            echo --input_file $path
        done) \
        $(if [[ {args[allow_seq_dict_incompatibility]} == "True" ]]; then
            echo --disable-sequence-dictionary-validation true
        fi)
}}
export -f run-shard

{run_shards}

{gather_shards}
"""
)
