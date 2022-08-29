# -*- coding: utf-8 -*-
# isort:skip_file
import os
import sys

from snakemake.shell import shell

# A hack is required for being able to import snappy_wrappers modules when in development mode.
# TODO: is there a more elegant way?
base_dir = os.path.normpath(os.path.join(os.path.dirname(__file__), "..", "..", "..", ".."))
sys.path.insert(0, base_dir)

from snappy_pipeline.utils import DictQuery


# Although optional for the tool, GATK recommend a providing a mappability track
map_bed = DictQuery(snakemake.config).get(
    "step_config/helper_gcnv_model/gcnv/path_uniquely_mapable_bed"
)


shell(
    r"""
set -x

gatk AnnotateIntervals \
    --interval-merging-rule OVERLAPPING_ONLY  \
    --mappability-track {map_bed} \
    --reference {snakemake.config[static_data_config][reference][path]} \
    --intervals {snakemake.input.interval_list} \
    --output {snakemake.output.tsv}
"""
)
