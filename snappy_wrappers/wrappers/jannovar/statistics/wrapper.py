from snappy_wrappers.snappy_wrapper import ShellWrapper

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

args = getattr(snakemake.params, "args", {})

ShellWrapper(snakemake).run(
    r"""
# See the following for the memory-related massaging
#
# http://bugs.java.com/view_bug.do?bug_id=8043516

# Call jannovar statistics
MALLOC_ARENA_MAX=4 \
jannovar \
    statistics \
    -XX:MaxHeapSize=4g \
    -XX:CompressedClassSpaceSize=1024m \
    -i {snakemake.input.vcf} \
    -o {snakemake.output.report} \
    -d {snakemake.input.path_ser}
"""
)
