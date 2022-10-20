# -*- coding: utf-8 -*-
# isort:skip_file
import os
import sys

from snakemake.shell import shell

# The following is required for being able to import snappy_wrappers modules
# inside wrappers.  These run in an "inner" snakemake process which uses its
# own conda environment which cannot see the snappy_pipeline installation.
base_dir = os.path.normpath(os.path.join(os.path.dirname(__file__), "..", "..", "..", ".."))
sys.path.insert(0, base_dir)

from snappy_pipeline.utils import DictQuery


# Although optional for the tool, GATK recommend a providing a mappability track
map_bed = DictQuery(snakemake.config).get(
    "step_config/helper_gcnv_model_wgs/gcnv/path_uniquely_mapable_bed"
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
