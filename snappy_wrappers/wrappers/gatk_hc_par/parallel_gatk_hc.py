# -*- coding: utf-8 -*-
"""Definition for GATK HC variant caller in parallel, genome is split into windows

isort:skip_file
"""

import os
import sys

# A hack is required for being able to import snappy_wrappers modules when in development mode.
base_dir = os.path.normpath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
sys.path.insert(0, base_dir)

from snappy_pipeline.workflows.abstract import ResourceUsage
from snappy_wrappers.wrapper_parallel import (
    ParallelVariantCallingBaseWrapper,
    gib,
    hours,
)


class ParallelGatkHcWrapper(ParallelVariantCallingBaseWrapper):
    """Parallel execution of GATK HaplotypeCaller

    Only the resource requirements are currently duplicated between the parallel GATK UG and
    GATK HC wrapper and these might need tweaking later on.  Thus, no shared parent class for
    GATK variant callers has been introduced (yet).
    """

    inner_wrapper = "gatk_hc"
    step_name = "variant_calling"
    tool_name = "gatk_hc"

    def __init__(self, snakemake):
        super().__init__(snakemake)
        self.job_resources = ResourceUsage(
            threads=1,
            memory=gib(14.0 * self.get_job_mult_memory()),
            time=str(hours(4 * self.get_job_mult_time())),
            partition=os.getenv('SNAPPY_PIPELINE_PARTITION'),
        )
        self.merge_resources = ResourceUsage(
            threads=2,
            memory=gib(2.0 * self.get_merge_mult_memory()),
            time=str(hours(4 * self.get_merge_mult_time())),
            partition=os.getenv('SNAPPY_PIPELINE_PARTITION'),
        )
