# -*- coding: utf-8 -*-
"""Wrapper for running Varscan for joint somatic/normal calling in parallel, genome is split into
windows
"""

import os
import sys
import textwrap

from snakemake import shell

# A hack is required for being able to import snappy_wrappers modules when in development mode.
# TODO: is there a more elegant way?
base_dir = os.path.normpath(os.path.join(os.path.dirname(__file__), "..", "..", "..", ".."))
sys.path.insert(0, base_dir)

from snappy_wrappers.wrapper_parallel import (
    gib,
    hours,
    ResourceUsage,
    ParallelVariantCallingBaseWrapper,
)


class ParallelVarscanCallJointWrapper(ParallelVariantCallingBaseWrapper):
    """Parallel execution of Varscan for somatic variant calling"""

    inner_wrapper = "varscan/call_joint"
    step_name = "somatic_variant_calling"
    tool_name = "varscan_joint"

    def __init__(self, snakemake):
        super().__init__(snakemake)
        self.job_resources = ResourceUsage(
            cores=2,
            memory=gib(3.75 * self.get_job_mult_memory()),
            duration=hours(4 * self.get_job_mult_time()),
        )
        self.merge_resources = ResourceUsage(
            cores=2,
            memory=gib(2.0 * self.get_merge_mult_memory()),
            duration=hours(4 * self.get_merge_mult_time()),
        )


# Write out information about conda installation.
shell(
    r"""
conda list >{snakemake.log.conda_list}
conda info >{snakemake.log.conda_info}
"""
)

# Kick off execution using the wrapper class defined above.
ParallelVarscanCallJointWrapper(snakemake).run()
