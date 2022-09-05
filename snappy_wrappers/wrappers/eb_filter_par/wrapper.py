# -*- coding: utf-8 -*-
"""Wrapper for running EasyBayes-Filter variant caller in parallel"""

from parallel_eb_filter import ParallelEasyBayesFilterWrapper
from snakemake import shell

# Kick off execution using the wrapper class defined above.
ParallelEasyBayesFilterWrapper(snakemake).run().shutdown_logging()

# Compute MD5 sums of logs.
shell(
    r"""
md5sum {snakemake.log.log} >{snakemake.log.log_md5}
"""
)
