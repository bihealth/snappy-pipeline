# -*- coding: utf-8 -*-
"""Wrapper for running Varscan for joint somatic/normal calling in parallel"""

from snakemake import shell

from parallel_call_joint import ParallelVarscanCallJointWrapper

# Write out information about conda installation.
shell(
    r"""
conda list >{snakemake.log.conda_list}
conda info >{snakemake.log.conda_info}
"""
)

# Kick off execution using the wrapper class defined above.
ParallelVarscanCallJointWrapper(snakemake).run()
