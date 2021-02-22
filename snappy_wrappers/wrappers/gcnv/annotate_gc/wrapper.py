# -*- coding: utf-8 -*-

import os
from snakemake.shell import shell

# A hack is required for being able to import snappy_wrappers modules when in development mode.
# TODO: is there a more elegant way?
base_dir = os.path.normpath(os.path.join(os.path.dirname(__file__), "..", "..", "..", ".."))
sys.path.insert(0, base_dir)

from snappy_pipeline.utils import DictQuery

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

map_bed = DictQuery(snakemake.config).get(
    "step_config/targeted_seq_cnv_calling/gcnv/path_uniquely_mapable_bed"
)

shell(
    r"""
set -x

gatk AnnotateIntervals \
    --interval-merging-rule OVERLAPPING_ONLY \
    --mappability-track {map_bed} \
    -R {snakemake.config[static_data_config][reference][path]} \
    -L {snakemake.input.interval_list} \
    -O {snakemake.output.tsv}
"""
)
