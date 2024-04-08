# -*- coding: utf-8 -*-
"""Wrapper for running Mutect2 variant caller in parallel"""

from parallel_mutect2 import ParallelMutect2Wrapper

# Kick off execution using the wrapper class defined above.
ParallelMutect2Wrapper(snakemake).run().shutdown_logging()
