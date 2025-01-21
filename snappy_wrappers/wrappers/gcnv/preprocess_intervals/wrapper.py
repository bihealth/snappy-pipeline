# -*- coding: utf-8 -*-
# isort:skip_file
from snappy_pipeline.utils import DictQuery
import os

from snakemake.shell import shell

args = getattr(snakemake.params, "args", {})

# Pick the target BED file to use.
for item in args["path_target_interval_list_mapping"]:
    if item["name"] == snakemake.wildcards.library_kit:
        target_interval_bed = item["path"]
        break
else:  # of for, did not break out
    raise Exception("Found no target intervals for %s" % snakemake.wildcards.library_kit)

shell(
    r"""
set -x

gatk PreprocessIntervals \
    --bin-length 0 \
    --interval-merging-rule OVERLAPPING_ONLY \
    -R {args[reference]} \
    -L {target_interval_bed} \
    -O {snakemake.output.interval_list}
"""
)
