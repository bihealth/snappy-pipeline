# -*- coding: utf-8 -*-
"""Execute Mutect variant caller with parallel wrapper."""

from snakemake import shell

from parallel_mutect_wrapper import ParallelMutectWrapper

# Kick off execution using the wrapper
ParallelMutectWrapper(snakemake).run().shutdown_logging()

# Compute MD5 sums of logs.
shell(
    r"""
md5sum {snakemake.log.log} >{snakemake.log.log_md5}
"""
)
