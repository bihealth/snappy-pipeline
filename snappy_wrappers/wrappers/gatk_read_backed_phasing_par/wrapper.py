# -*- coding: utf-8 -*-
"""Wrapper code for GATK ReadBackedPhasing

isort:skip_file
"""

import os
import sys
import textwrap

from snakemake.shell import shell

# A hack is required for being able to import snappy_wrappers modules when in development mode.
# TODO: is there a more elegant way?
base_dir = os.path.normpath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
sys.path.insert(0, base_dir)

from snappy_wrappers.wrapper_parallel import (
    ParallelVariantAnnotationBaseWrapper,
    ResourceUsage,
    gib,
    hours,
)

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>"


class ParallelGaktReadBackedPhasingWrapper(ParallelVariantAnnotationBaseWrapper):
    """Parallel execution of GATK ReadBackedPhasings"""

    inner_wrapper = "gatk_read_backed_phasing"
    step_name = "variant_phasing"
    tool_name = "gatk_read_backed_phasing"
    forward_input_keys = ("vcf", "tbi", "bam")

    def __init__(self, snakemake):
        super().__init__(snakemake)
        self.job_resources = ResourceUsage(
            cores=2,
            memory=gib(14.0 * self.get_job_mult_memory()),
            duration=hours(4 * self.get_job_mult_time()),
        )
        self.merge_resources = ResourceUsage(
            cores=2,
            memory=gib(2.0 * self.get_merge_mult_memory()),
            duration=hours(4 * self.get_merge_mult_time()),
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
                        for key in ("vcf", "tbi", "bam")
                    }
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
                        **{input_},
                    output:
                        touch("job_out.{jobno}.d/.done"),
                        **{output}
                    params:
                        **{params}
                    wrapper: '{wrapper_prefix}/snappy_wrappers/wrappers/{inner_wrapper}'

                cluster_config['chunk_{jobno}'] = {resources}
            """
            ).format(**vals).lstrip()


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
ParallelGaktReadBackedPhasingWrapper(snakemake).run()

# Compute MD5 sums of logs.
shell(
    r"""
md5sum {snakemake.log.log} >{snakemake.log.log_md5}
"""
)
