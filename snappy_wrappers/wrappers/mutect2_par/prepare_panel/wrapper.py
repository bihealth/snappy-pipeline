# -*- coding: utf-8 -*-
"""Wrapper for running MuTect 2 variant caller in parallel, genome is split into windows

isort:skip_file
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
    ParallelVariantCallingBaseWrapper,
    ResourceUsage,
    gib,
    hours,
)


class ParallelMutect2Wrapper(ParallelVariantCallingBaseWrapper):
    """Parallel execution of MuTect 2"""

    inner_wrapper = "mutect2/prepare_panel"
    step_name = "panel_of_normals"
    tool_name = "mutect2"

    def __init__(self, snakemake):
        super().__init__(snakemake)
        self.job_resources = ResourceUsage(
            cores=1,
            memory=gib(14.0 * self.get_job_mult_memory()),
            duration=hours(4 * self.get_job_mult_time()),
        )
        self.merge_resources = ResourceUsage(
            cores=1,
            memory=gib(2.0 * self.get_merge_mult_memory()),
            duration=hours(4 * self.get_merge_mult_time()),
        )

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
                "input_bam": repr(self.snakemake.input.normal_bam),
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
                        {input_bam},
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
ParallelMutect2Wrapper(snakemake).run().shutdown_logging()

# Compute MD5 sums of logs.
shell(
    r"""
md5sum {snakemake.log.log} >{snakemake.log.log_md5}
"""
)
