# -*- coding: utf-8 -*-
"""Wrapper for running MuTect 2 variant caller in parallel, genome is split into windows

isort:skip_file
"""

import os
import sys

from snakemake import shell

from parallel_mutect2_wrapper import ParallelMutect2Wrapper

# A hack is required for being able to import snappy_wrappers modules when in development mode.
# TODO: is there a more elegant way?
base_dir = os.path.normpath(os.path.join(os.path.dirname(__file__), "..", "..", "..", ".."))
sys.path.insert(0, base_dir)


# Kick off execution using the wrapper class defined above.
ParallelMutect2Wrapper(snakemake).run().shutdown_logging()

# Compute MD5 sums of logs.
shell(
    r"""
md5sum {snakemake.log.log} >{snakemake.log.log_md5}
"""
)
