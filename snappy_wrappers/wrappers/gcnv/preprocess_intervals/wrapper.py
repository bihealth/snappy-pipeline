# -*- coding: utf-8 -*-
# isort:skip_file

from typing import TYPE_CHECKING

from snappy_wrappers.snappy_wrapper import ShellWrapper

if TYPE_CHECKING:
    from snakemake.iocontainers import snakemake

# Targeted sequencing: restrict the bins to the target regions of the library kit
target_bed = snakemake.input.get("target_bed", "")
target_interval_bed = f"-L {target_bed}" if target_bed else ""

ShellWrapper(snakemake).run(
    r"""

gatk PreprocessIntervals \
    --bin-length 0 \
    --interval-merging-rule OVERLAPPING_ONLY \
    -R {snakemake.input.reference} \
    {target_interval_bed} \
    -O {snakemake.output.interval_list}
"""
)
