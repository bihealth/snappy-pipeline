# -*- coding: utf-8 -*-
"""Definition for running GATK UG variant caller in parallel, genome is split into windows"""

import os
import sys

# The following is required for being able to import snappy_wrappers modules
# inside wrappers.  These run in an "inner" snakemake process which uses its
# own conda environment which cannot see the snappy_pipeline installation.
base_dir = os.path.normpath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
sys.path.insert(0, base_dir)

from snappy_wrappers.resource_usage import ResourceUsage  # noqa: E402
from snappy_wrappers.wrapper_parallel import (  # noqa: E402
    ParallelVariantCallingBaseWrapper,
    gib_to_string,
    hours,
)


class ParallelGatkUgWrapper(ParallelVariantCallingBaseWrapper):
    """Parallel execution of GATK UnifiedGenotyper

    Only the resource requirements are currently duplicated between the parallel GATK UG and
    GATK HC wrapper and these might need tweaking later on. Thus, no shared parent class for
    GATK variant callers has been introduced (yet).
    """

    inner_wrapper = "gatk_ug"
    step_name = "variant_calling"
    tool_name = "gatk_ug"

    def __init__(self, snakemake):
        super().__init__(snakemake)
        self.job_resources = ResourceUsage(
            threads=1,
            memory=gib_to_string(14.0 * self.get_job_mult_memory()),
            time=hours(4 * self.get_job_mult_time()),
            partition=os.getenv("SNAPPY_PIPELINE_PARTITION"),
        )
        self.merge_resources = ResourceUsage(
            threads=2,
            memory=gib_to_string(2.0 * self.get_merge_mult_memory()),
            time=hours(4 * self.get_merge_mult_time()),
            partition=os.getenv("SNAPPY_PIPELINE_PARTITION"),
        )
