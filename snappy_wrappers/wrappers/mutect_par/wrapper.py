# -*- coding: utf-8 -*-
"""Wrapper for running Mutect2variant caller in parallel"""

from parallel_mutect import ParallelMutectWrapper
from snakemake import shell

# Kick off execution using the wrapper class defined above.
ParallelMutectWrapper(snakemake).run().shutdown_logging()

# Compute MD5 sums of logs.
shell(
    r"""
md5sum {snakemake.log.log} >{snakemake.log.log_md5}
"""
)
