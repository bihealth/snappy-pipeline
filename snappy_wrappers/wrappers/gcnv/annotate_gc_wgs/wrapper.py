# -*- coding: utf-8 -*-
# isort:skip_file
from snakemake.shell import shell

args = getattr(snakemake.params, "args", {})

# Although optional for the tool, GATK recommend a providing a mappability track
map_bed = args["path_uniquely_mapable_bed"]

shell(
    r"""
set -x

gatk AnnotateIntervals \
    --interval-merging-rule OVERLAPPING_ONLY  \
    --mappability-track {map_bed} \
    --reference {args[reference]} \
    --intervals {snakemake.input.interval_list} \
    --output {snakemake.output.tsv}
"""
)
