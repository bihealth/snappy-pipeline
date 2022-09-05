# -*- coding: utf-8 -*-
"""Wrapper for Jannovar germline annotation in parallel"""

from parallel_annotate_vcf import ParallelJannovarAnnotateVcfWrapper
from snakemake import shell

# Kick off execution using the wrapper class defined above.
ParallelJannovarAnnotateVcfWrapper(snakemake).run().shutdown_logging()

# Compute MD5 sums of logs.
shell(
    r"""
md5sum {snakemake.log.log} >{snakemake.log.log_md5}
"""
)
