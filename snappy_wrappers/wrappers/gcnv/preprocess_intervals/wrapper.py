# -*- coding: utf-8 -*-
# isort:skip_file
from snappy_pipeline.utils import DictQuery
import os

from snakemake.shell import shell


# Pick the target BED file to use.
config = DictQuery(snakemake.config).get("step_config/sv_calling_targeted/gcnv")
for item in config["path_target_interval_list_mapping"]:
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
    -R {snakemake.config[static_data_config][reference][path]} \
    -L {target_interval_bed} \
    -O {snakemake.output.interval_list}
"""
)
