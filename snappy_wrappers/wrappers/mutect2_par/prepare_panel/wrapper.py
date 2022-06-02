# -*- coding: utf-8 -*-
"""Wrapper for running Mutect2 variant caller in parallel"""

from snakemake import shell

from parallel_prepare_panel import ParallelMutect2Wrapper

# Kick off execution using the wrapper class defined above.
ParallelMutect2Wrapper(snakemake).run().shutdown_logging()

# Compute MD5 sums of logs.
shell(
    r"""
md5sum {snakemake.log.log} >{snakemake.log.log_md5}
"""
)
