# -*- coding: utf-8 -*-
"""Definition for GATK ReadBackedPhasing paralllel wrapper
"""

import os
import sys
import textwrap

# The following is required for being able to import snappy_wrappers modules
# inside wrappers.  These run in an "inner" snakemake process which uses its
# own conda environment which cannot see the snappy_pipeline installation.
base_dir = os.path.normpath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
sys.path.insert(0, base_dir)

from snappy_wrappers.resource_usage import ResourceUsage  # noqa: E402
from snappy_wrappers.wrapper_parallel import (  # noqa: E402
    ParallelVariantAnnotationBaseWrapper,
    gib_to_string,
    hours,
)

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"


class ParallelGaktReadBackedPhasingWrapper(ParallelVariantAnnotationBaseWrapper):
    """Parallel execution of GATK ReadBackedPhasings"""

    inner_wrapper = "gatk_read_backed_phasing"
    step_name = "variant_phasing"
    tool_name = "gatk_read_backed_phasing"
    forward_input_keys = ("vcf", "vcf_tbi", "bam")

    def __init__(self, snakemake):
        super().__init__(snakemake)
        self.job_resources = ResourceUsage(
            threads=2,
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

    # TODO: this is a clever trick and could go into super class
    def _abs_path(self, path):
        """Helper to build absolute paths"""
        if isinstance(path, str):
            return os.path.realpath(os.path.join(self.main_cwd, path))
        else:
            return list(map(self._abs_path, list(path)))

    def construct_parallel_rules(self):
        """Construct the rules for parallel processing to generate."""
        for jobno, region in enumerate(self.get_regions()):
            params = dict(self.snakemake.params)
            params.setdefault("args", {}).update({"intervals": [region.human_readable()]})
            output = {
                key: "job_out.{jobno}.d/out/tmp_{jobno}.{ext}".format(jobno=jobno, ext=ext)
                for key, ext in self.key_ext.items()
            }
            vals = {
                "input_": repr(
                    {
                        key: self._abs_path(getattr(self.snakemake.input, key))
                        for key in ("vcf", "vcf_tbi", "bam")
                    }
                ),
                "jobno": jobno,
                "params": repr(params),
                "output": repr(output),
                "wrapper_prefix": "file://" + self.wrapper_base_dir,
                "inner_wrapper": self.inner_wrapper,
            }
            yield textwrap.dedent(
                r"""
                rule chunk_{jobno}:
                    input:
                        **{input_},
                    output:
                        touch("job_out.{jobno}.d/.done"),
                        **{output}
                    threads: resource_merge_threads
                    resources:
                        time=resource_merge_time,
                        memory=resource_merge_memory,
                        partition=resource_merge_partition,
                    params:
                        **{params}
                    wrapper: '{wrapper_prefix}/snappy_wrappers/wrappers/{inner_wrapper}'

            """
            ).format(**vals).lstrip()
