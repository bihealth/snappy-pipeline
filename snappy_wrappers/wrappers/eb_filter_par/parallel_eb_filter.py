# -*- coding: utf-8 -*-
"""Definition for running EasyBayes-Filter in parallel, genome is split into windows"""

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
    ParallelVcfOutputBaseWrapper,
    gib_to_string,
    hours,
)


class ParallelEasyBayesFilterWrapper(ParallelVcfOutputBaseWrapper):
    """Parallel execution of somatic variant annotation: ``eb_filter``."""

    inner_wrapper = "eb_filter"
    step_name = "somatic_variant_filtration"
    tool_name = "eb_filter"

    def __init__(self, snakemake):
        super().__init__(snakemake)
        self.job_resources = ResourceUsage(
            threads=2,
            memory=gib_to_string(8.0 * self.get_job_mult_memory()),
            time=hours(4 * self.get_job_mult_time()),
            partition=os.getenv("SNAPPY_PIPELINE_PARTITION"),
        )
        self.merge_resources = ResourceUsage(
            threads=2,
            memory=gib_to_string(2.0 * self.get_merge_mult_memory()),
            time=hours(4 * self.get_merge_mult_time()),
            partition=os.getenv("SNAPPY_PIPELINE_PARTITION"),
        )

    def construct_parallel_rules(self):
        """Construct the rules for parallel processing to generate."""
        for jobno, region in enumerate(self.get_regions()):
            params = {}
            for p in self.snakemake.params:
                params.update(p)
            params.setdefault("args", {}).update({"interval": region.human_readable(False)})
            output = {
                key: "job_out.{jobno}.d/out/tmp_{jobno}.{ext}".format(jobno=jobno, ext=ext)
                for key, ext in self.key_ext.items()
            }
            vals = {
                "input_vcf": repr(
                    os.path.realpath(os.path.join(self.main_cwd, self.snakemake.input.vcf))
                ),
                "input_bam": repr(
                    os.path.realpath(os.path.join(self.main_cwd, self.snakemake.input.bam))
                ),
                "input_txt": repr(
                    os.path.realpath(os.path.join(self.main_cwd, self.snakemake.input.txt))
                ),
                "jobno": jobno,
                "params": repr(params),
                "output": repr(output),
                "wrapper_prefix": "file://" + self.wrapper_base_dir,
                "inner_wrapper": self.inner_wrapper,
            }
            yield (
                textwrap.dedent(
                    r"""
                rule chunk_{jobno}:
                    input:
                        vcf={input_vcf},
                        bam={input_bam},
                        txt={input_txt},
                    output:
                        touch("job_out.{jobno}.d/.done"),
                        **{output}
                    threads: resource_chunk_threads
                    resources:
                        time=resource_chunk_time,
                        memory=resource_chunk_memory,
                        partition=resource_chunk_partition,
                    params:
                        **{params}
                    wrapper: '{wrapper_prefix}/snappy_wrappers/wrappers/{inner_wrapper}'


            """
                )
                .format(**vals)
                .lstrip()
            )
