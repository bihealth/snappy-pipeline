# -*- coding: utf-8 -*-
# isort:skip_file
import os
import sys

from snakemake.shell import shell

# Although optional for the tool, GATK recommend a providing a mappability track
map_bed = snakemake.config["step_config"]["helper_gcnv_model_wgs"]["gcnv"]["path_uniquely_mapable_bed"]
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
