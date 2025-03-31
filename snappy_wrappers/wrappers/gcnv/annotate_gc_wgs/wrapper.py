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

args = getattr(snakemake.params, "args", {})

from snappy_pipeline.utils import DictQuery


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
