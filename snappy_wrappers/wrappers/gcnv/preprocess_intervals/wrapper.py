# -*- coding: utf-8 -*-
# isort:skip_file
from snappy_pipeline.utils import DictQuery
import os

from snakemake.shell import shell


# Pick the target BED file to use.  If it goes by the name "default" then we
# simply take the one at xhmm/path_target_interval_list.  Otherwise, we
# have to go through the list in xhmm/path_target_interval_list_mapping.
gcnv_config = DictQuery(snakemake.config).get("step_config/targeted_seq_cnv_calling/gcnv")
if snakemake.wildcards.library_kit == "default":
    target_interval_bed = gcnv_config["path_target_interval_list"]
else:
    for item in gcnv_config["path_target_interval_list_mapping"]:
        if item["name"] == snakemake.wildcards.library_kit:
            target_interval_bed = item["path"]
            break
    else:  # of for, did not break out
        raise Exception("Found no target intervals for %s" % item["name"])

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
