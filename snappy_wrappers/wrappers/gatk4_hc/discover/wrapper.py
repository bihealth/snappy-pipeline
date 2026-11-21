from snappy_wrappers.snappy_wrapper import ShellWrapper
from snappy_wrappers.tools import shards
from snappy_wrappers.tools.java import java_heap_mb

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

args = getattr(snakemake.params, "args", {})

# Each of the parallel shards gets an equal part of 75% of the job memory as Java heap
java_mem_mb = java_heap_mb(snakemake, shards=args["num_threads"])
annotation_groups = " ".join(f"-G {group}" for group in args["annotation_groups"])
allow_seq_dict_incompatibility = (
    "--disable-sequence-dictionary-validation true"
    if args.get("allow_seq_dict_incompatibility", False)
    else ""
)
shard_intervals = shards.intervals(
    snakemake.input.reference,
    snakemake.input.reference_fai,
    args["window_length"],
    args["padding"],
    args["ignore_contigs"],
)
run_shards = shards.run_shards(args["num_threads"])
gather_shards = shards.gather(snakemake.input.reference, snakemake.output.gvcf)

ShellWrapper(snakemake).run(
    r"""
{shard_intervals}

# Function to run HaplotypeCaller in to generate gVCF
run-shard() {{
    job_no=$1
    interval=$2

    # Create per-sample gVCF file
    GATK_JAVA_MEMORY={java_mem_mb}m
    gatk \
        HaplotypeCaller \
        --java-options "-Xmx$GATK_JAVA_MEMORY -Djava.io.tmpdir=$TMPDIR" \
        --tmp-dir $TMPDIR \
        --output $TMPDIR/shards-output/$(printf %06d $job_no).g.vcf.gz \
        --reference {snakemake.input.reference} \
        --dbsnp {snakemake.input.dbsnp} \
        --intervals $interval \
        --input {snakemake.input.bam} \
        {annotation_groups} \
        {allow_seq_dict_incompatibility} \
        {extra_arguments} \
        -ERC GVCF
}}
export -f run-shard

{run_shards}

{gather_shards}
"""
)
