# -*- coding: utf-8 -*-
# isort:skip_file
from snappy_pipeline.utils import DictQuery

from snakemake.shell import shell

# Pick the target BED file to use.
# FIXME: why is 'target_interval_bed' not used?
config = DictQuery(snakemake.config).get("step_config/helper_gcnv_model_targeted/gcnv")
for item in config["path_target_interval_list_mapping"]:
    if item["name"] == snakemake.wildcards.library_kit:
        target_interval_bed = item["path"]
        break
else:  # of for, did not break out
    raise Exception("Found no target intervals for %s" % item["name"])

map_bed = DictQuery(snakemake.config).get(
    "step_config/helper_gcnv_model_targeted/gcnv/path_uniquely_mapable_bed"
)

shell(
    r"""
set -x

gatk AnnotateIntervals \
    --interval-merging-rule OVERLAPPING_ONLY \
    --mappability-track {map_bed} \
    --reference {snakemake.config[static_data_config][reference][path]} \
    --intervals {snakemake.input.interval_list} \
    --output {snakemake.output.tsv}
"""
)
