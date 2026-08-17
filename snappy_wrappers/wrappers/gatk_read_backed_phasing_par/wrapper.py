# -*- coding: utf-8 -*-
"""Wrapper code for GATK ReadBackedPhasing"""

from typing import TYPE_CHECKING

from parallel_read_backed_phasing import ParallelGaktReadBackedPhasingWrapper
from snappy_wrappers.snappy_wrapper import ShellWrapper

if TYPE_CHECKING:
    from snakemake.iocontainers import snakemake

# Kick off execution using the wrapper class defined above.
ParallelGaktReadBackedPhasingWrapper(snakemake).run()

# Trigger standardized conda/log capture and md5 generation from SnappyWrapper.
ShellWrapper(snakemake, with_output_links=False).run("true")
