# -*- coding: utf-8 -*-
"""Definition for Jannovar somatic annotation in parallel, genome is split into windows"""

import os
import sys

# The following is required for being able to import snappy_wrappers modules
# inside wrappers.  These run in an "inner" snakemake process which uses its
# own conda environment which cannot see the snappy_pipeline installation.
base_dir = os.path.normpath(os.path.join(os.path.dirname(__file__), "..", "..", "..", ".."))
sys.path.insert(0, base_dir)

from snappy_wrappers.resource_usage import ResourceUsage  # noqa: E402
from snappy_wrappers.wrapper_parallel import (  # noqa: E402
    ParallelSomaticVariantAnnotationBaseWrapper,
    gib_to_string,
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
            threads=2,
            memory=gib_to_string(10.0 * self.get_job_mult_memory()),
            time=hours(4 * self.get_job_mult_time()),
            partition=os.getenv("SNAPPY_PIPELINE_PARTITION"),
        )
        self.merge_resources = ResourceUsage(
            threads=2,
            memory=gib_to_string(2.0 * self.get_merge_mult_memory()),
            time=hours(4 * self.get_merge_mult_time()),
            partition=os.getenv("SNAPPY_PIPELINE_PARTITION"),
        )
