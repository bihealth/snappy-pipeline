# -*- coding: utf-8 -*-
# isort:skip_file
from typing import TYPE_CHECKING

from snappy_wrappers.snappy_wrapper import ShellWrapper

if TYPE_CHECKING:
    from snakemake.iocontainers import snakemake

# Although optional for the tool, GATK recommend a providing a mappability track
ShellWrapper(snakemake).run(
    r"""

gatk AnnotateIntervals \
    --interval-merging-rule OVERLAPPING_ONLY  \
    --mappability-track {snakemake.input.mappability} \
    --reference {snakemake.input.reference} \
    --intervals {snakemake.input.interval_list} \
    --output {snakemake.output.tsv}
"""
)
