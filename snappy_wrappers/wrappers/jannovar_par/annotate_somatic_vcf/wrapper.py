# -*- coding: utf-8 -*-
"""Wrapper for Jannovar somatic annotation in parallel"""

from snakemake import shell

from parallel_annotate_somatic_vcf import ParallelJannovarAnnotateSomaticVcfWrapper

# Kick off execution using the wrapper class defined above.
ParallelJannovarAnnotateSomaticVcfWrapper(snakemake).run().shutdown_logging()

# Compute MD5 sums of logs.
shell(
    r"""
md5sum {snakemake.log.log} >{snakemake.log.log_md5}
"""
)
