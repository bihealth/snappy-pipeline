# -*- coding: utf-8 -*-
"""Definition for Jannovar somatic annotation in parallel, genome is split into windows

isort:skip_file
"""

from snappy_wrappers.wrapper_parallel import (
    ParallelSomaticVariantAnnotationBaseWrapper,
    gib_to_string,
    hours,
)
from snappy_wrappers.resource_usage import ResourceUsage
import os


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
