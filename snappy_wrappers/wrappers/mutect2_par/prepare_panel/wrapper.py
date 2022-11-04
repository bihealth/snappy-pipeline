# -*- coding: utf-8 -*-
"""Wrapper for running Mutect2 variant caller in parallel"""

import os
import sys

# The following is required for being able to import snappy_wrappers modules
# inside wrappers.  These run in an "inner" snakemake process which uses its
# own conda environment which cannot see the snappy_pipeline installation.
base_dir = os.path.normpath(os.path.join(os.path.dirname(__file__), "..", "..", "..", ".."))
sys.path.insert(0, base_dir)

from parallel_prepare_panel import ParallelMutect2Wrapper
from snakemake import shell

# Kick off execution using the wrapper class defined above.
ParallelMutect2Wrapper(snakemake).run().shutdown_logging()

# Compute MD5 sums of logs.
shell(
    r"""
md5sum {snakemake.log.log} >{snakemake.log.log_md5}
"""
)
