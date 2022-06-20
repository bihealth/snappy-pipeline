# -*- coding: utf-8 -*-
"""Wrapper code for GATK ReadBackedPhasing"""

from snakemake.shell import shell

from parallel_read_backed_phasing import ParallelGaktReadBackedPhasingWrapper

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
ParallelGaktReadBackedPhasingWrapper(snakemake).run()

# Compute MD5 sums of logs.
shell(
    r"""
md5sum {snakemake.log.log} >{snakemake.log.log_md5}
"""
)
