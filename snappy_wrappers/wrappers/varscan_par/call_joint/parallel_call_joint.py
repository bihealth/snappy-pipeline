# -*- coding: utf-8 -*-
"""Definition for running Varscan for joint somatic/normal calling in parallel, genome is split into
windows

isort:skip_file
"""

import os


from snappy_wrappers.resource_usage import ResourceUsage  # noqa: E402
from snappy_wrappers.wrapper_parallel import (  # noqa: E402
    ParallelVariantCallingBaseWrapper,
    gib_to_string,
    hours,
)


class ParallelVarscanCallJointWrapper(ParallelVariantCallingBaseWrapper):
    """Parallel execution of Varscan for somatic variant calling"""

    inner_wrapper = "varscan/call_joint"
    step_name = "somatic_variant_calling"
    tool_name = "varscan_joint"

    def __init__(self, snakemake):
        super().__init__(snakemake)
        self.job_resources = ResourceUsage(
            threads=2,
            memory=gib_to_string(3.75 * self.get_job_mult_memory()),
            time=hours(4 * self.get_job_mult_time()),
            partition=os.getenv("SNAPPY_PIPELINE_PARTITION"),
        )
        self.merge_resources = ResourceUsage(
            threads=2,
            memory=gib_to_string(2.0 * self.get_merge_mult_memory()),
            time=hours(4 * self.get_merge_mult_time()),
            partition=os.getenv("SNAPPY_PIPELINE_PARTITION"),
        )
