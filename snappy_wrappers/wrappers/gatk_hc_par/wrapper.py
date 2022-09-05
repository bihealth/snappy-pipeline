# -*- coding: utf-8 -*-
"""Wrapper for running GATK HC variant caller in parallel
"""

import os

from parallel_gatk_hc import ParallelGatkHcWrapper
from snakemake import io, shell

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
snakemake.input = io.Namedlist(map(os.path.realpath, snakemake.input))
ParallelGatkHcWrapper(snakemake).run().shutdown_logging()

# Compute MD5 sums of logs.
shell(
    r"""
md5sum {snakemake.log.log} >{snakemake.log.log_md5}
"""
)
