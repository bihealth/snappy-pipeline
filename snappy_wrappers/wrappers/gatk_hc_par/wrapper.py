# -*- coding: utf-8 -*-
"""Wrapper for running GATK HC variant caller in parallel, genome is split into windows

isort:skip_file
"""

import os
import sys

from snakemake import io, shell

# A hack is required for being able to import snappy_wrappers modules when in development mode.
# TODO: is there a more elegant way?
base_dir = os.path.normpath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
sys.path.insert(0, base_dir)

from snappy_wrappers.wrapper_parallel import (
    ParallelVariantCallingBaseWrapper,
    ResourceUsage,
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
    step_name = snakemake.params.step_key
    tool_name = snakemake.params.caller_key

    def __init__(self, snakemake):
        super().__init__(snakemake)
        self.job_resources = ResourceUsage(
            cores=1,
            memory=gib(14.0 * self.get_job_mult_memory()),
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
md5sum {snakemake.log.conda_list} >{snakemake.log.conda_list_md5}
md5sum {snakemake.log.conda_info} >{snakemake.log.conda_info_md5}
"""
)

# Kick off execution using the wrapper class defined above.
snakemake.input = io.Namedlist(map(os.path.realpath, snakemake.input))
ParallelGatkHcWrapper(snakemake).run().shutdown_logging()

# Compute MD5 sums of logs.
shell(
    r"""
md5sum {snakemake.log.log} >{snakemake.log.log_md5}
"""
)
