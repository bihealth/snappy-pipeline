# -*- coding: utf-8 -*-
"""Execute Mutect2 variant caller with parallel wrapper."""

from snakemake import shell

from parallel_mutect2_wrapper import ParallelMutect2Wrapper

# Kick off execution using the wrapper class
ParallelMutect2Wrapper(snakemake).run().shutdown_logging()

# Compute MD5 sums of logs.
shell(
    r"""
md5sum {snakemake.log.log} >{snakemake.log.log_md5}
"""
)
