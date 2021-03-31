# -*- coding: utf-8 -*-
"""Wrapper for Jannovar somatic annotation in parallel, genome is split into windows

isort:skip_file
"""

import os
import sys

from snakemake import shell

# A hack is required for being able to import snappy_wrappers modules when in development mode.
# TODO: is there a more elegant way?
base_dir = os.path.normpath(os.path.join(os.path.dirname(__file__), "..", "..", "..", ".."))
sys.path.insert(0, base_dir)

from snappy_wrappers.wrapper_parallel import (
    ParallelSomaticVariantAnnotationBaseWrapper,
    ResourceUsage,
    gib,
    hours,
)


class ParallelJannovarAnnotateSomaticVcfWrapper(ParallelSomaticVariantAnnotationBaseWrapper):
    """Parallel execution of somatic ``jannovar annotate-vcf``."""

    inner_wrapper = "jannovar/annotate_somatic_vcf"
    step_name = "somatic_variant_annotation"
    tool_name = None  # tool confg is on top level

    def __init__(self, snakemake):
        super().__init__(snakemake)
        self.job_resources = ResourceUsage(
            cores=2,
            memory=gib(10.0 * self.get_job_mult_memory()),
            duration=hours(4 * self.get_job_mult_time()),
        )
        self.merge_resources = ResourceUsage(
            cores=2,
            memory=gib(2.0 * self.get_merge_mult_memory()),
            duration=hours(4 * self.get_merge_mult_time()),
        )


# Kick off execution using the wrapper class defined above.
ParallelJannovarAnnotateSomaticVcfWrapper(snakemake).run().shutdown_logging()

# Compute MD5 sums of logs.
shell(
    r"""
md5sum {snakemake.log.log} >{snakemake.log.log_md5}
"""
)
