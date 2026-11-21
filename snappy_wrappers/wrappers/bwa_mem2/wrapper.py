from typing import TYPE_CHECKING

from snappy_wrappers.snappy_wrapper import ShellWrapper
from snappy_wrappers.tools import bwa_mem

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

read_functions = bwa_mem.read_functions(
    input_left, input_right, read_groups, trim_adapters, num_threads_trimming, mask_duplicates
)
mapping_pipeline = bwa_mem.mapping_pipeline(
    "run_bwa_mem2",
    snakemake.output,
    num_threads_bam_view,
    memory_bam_sort,
    num_threads_bam_sort,
)

ShellWrapper(snakemake).run(
    r"""
{read_functions}

# Run BWA-MEM2 on the reads of read group $1
run_bwa_mem2() {{
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

{mapping_pipeline}
"""
)
