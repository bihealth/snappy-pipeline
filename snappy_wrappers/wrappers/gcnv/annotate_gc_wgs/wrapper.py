# -*- coding: utf-8 -*-
# isort:skip_file
from snakemake.shell import shell

# Although optional for the tool, GATK recommend a providing a mappability track
map_bed = snakemake.config['step_config'][snakemake.params.step_key]['gcnv']['path_uniquely_mapable_bed']
# Note: this wrapper is only used in ther model building workflow, not in calling from a built model

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
