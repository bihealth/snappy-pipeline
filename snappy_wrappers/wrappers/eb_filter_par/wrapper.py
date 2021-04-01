# -*- coding: utf-8 -*-
"""Wrapper for running EasyBayes-Filter in parallel, genome is split into windows

isort:skip_file
"""

import os
import sys
import textwrap

from snakemake import shell

# A hack is required for being able to import snappy_wrappers modules when in development mode.
# TODO: is there a more elegant way?
base_dir = os.path.normpath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
sys.path.insert(0, base_dir)

# pylint: disable=E402
from snappy_wrappers.wrapper_parallel import (
    ParallelVcfOutputBaseWrapper,
    ResourceUsage,
    gib,
    hours,
)
# pylint: enable=E402

class ParallelEasyBayesFilterWrapper(ParallelVcfOutputBaseWrapper):
    """Parallel execution of somatic variant annotation: eb_filter``."""

    inner_wrapper = "eb_filter"
    step_name = "somatic_variant_filtration"
    tool_name = "eb_filter"

    def __init__(self, snakemake):
        super().__init__(snakemake)
        self.job_resources = ResourceUsage(
            cores=2,
            memory=gib(8.0 * self.get_job_mult_memory()),
            duration=hours(4 * self.get_job_mult_time()),
        )
        self.merge_resources = ResourceUsage(
            cores=2,
            memory=gib(2.0 * self.get_merge_mult_memory()),
            duration=hours(4 * self.get_merge_mult_time()),
        )

    def construct_parallel_rules(self):
        """Construct the rules for parallel processing to generate."""
        for jobno, region in enumerate(self.get_regions()):
            params = dict(self.snakemake.params)
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
                "resources": repr(self.res_converter(self.job_resources).to_res_dict()),
            }
            yield textwrap.dedent(
                r"""
                rule chunk_{jobno}:
                    input:
                        vcf={input_vcf},
                        bam={input_bam},
                        txt={input_txt},
                    output:
                        touch("job_out.{jobno}.d/.done"),
                        **{output}
                    params:
                        **{params}
                    wrapper: '{wrapper_prefix}/snappy_wrappers/wrappers/{inner_wrapper}'


                cluster_config['chunk_{jobno}'] = {resources}
            """
            ).format(**vals).lstrip()


# Kick off execution using the wrapper class defined above.
ParallelEasyBayesFilterWrapper(snakemake).run().shutdown_logging()

# Compute MD5 sums of logs.
shell(
    r"""
md5sum {snakemake.log.log} >{snakemake.log.log_md5}
"""
)
