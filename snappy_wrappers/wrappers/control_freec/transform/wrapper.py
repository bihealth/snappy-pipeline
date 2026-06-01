# -*- coding: utf-8 -*-
"""Wrapper for merging multiple tables in R on shared columns"""

import os

from snakemake import shell

from snappy_wrappers.snappy_wrapper import ShellWrapper

shell.executable("/bin/bash")

args = getattr(snakemake.params, "args", {})

rscript = os.path.join(
    os.path.dirname(os.path.realpath(__file__)), "snappy-convert-control_freec.R"
)

ShellWrapper(snakemake).run(
    r"""
set -x

R --vanilla -e "source(\"{rscript}\") ; library(magrittr) ; \
    control_freec_write_files( \
    sample_name = \"{args[cancer_library]}\", \
    ratios_fn = \"{snakemake.input.ratio}\", \
    log2_fn = \"{snakemake.output.log2}\", \
    call_fn = \"{snakemake.output.call}\", \
    segments_fn = \"{snakemake.output.segments}\", \
    cns_fn = \"{snakemake.output.cns}\", \
    cnr_fn = \"{snakemake.output.cnr}\", \
    org_obj={args[org_obj]}, \
    tx_obj={args[tx_obj]}, \
    bs_obj={args[bs_obj]})"


"""
)
