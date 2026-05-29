# -*- coding: utf-8 -*-
# isort:skip_file
from typing import TYPE_CHECKING

from snappy_wrappers.snappy_wrapper import ShellWrapper

if TYPE_CHECKING:
    from snakemake.iocontainers import snakemake

args = getattr(snakemake.params, "args", {})

# Pick the target BED file to use.
# FIXME: why is 'target_interval_bed' not used?
for item in args["path_target_interval_list_mapping"]:
    if item["name"] == snakemake.wildcards.library_kit:
        target_interval_bed = item["path"]
        break
else:  # of for, did not break out
    raise Exception("Found no target intervals for %s" % item["name"])

map_bed = args["path_uniquely_mapable_bed"]

ShellWrapper(snakemake).run(
    r"""

gatk AnnotateIntervals \
    --interval-merging-rule OVERLAPPING_ONLY \
    --mappability-track {map_bed} \
    --reference {args[reference]} \
    --intervals {snakemake.input.interval_list} \
    --output {snakemake.output.tsv}
"""
)
