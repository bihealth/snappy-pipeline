# -*- coding: utf-8 -*-
"""Wrapper for merging multiple tables in R on shared columns"""

import os

from snappy_wrappers.snappy_wrapper import ShellWrapper

__author__ = "Eric Blanc <eric.blanc@bih-charite.de>"

r_script = os.path.abspath(os.path.join(os.path.dirname(__file__), "script.R"))
helper_functions = os.path.join(os.path.dirname(r_script), "..", "helper_functions.R")

args = getattr(snakemake.params, "args", {})

ShellWrapper(snakemake).run(
    r"""
# Run the R script --------------------------------------------------------------------------------

R --vanilla --slave << __EOF
source("{helper_functions}")
source("{r_script}")
write.table(
    cns_to_cna("{snakemake.input.DNAcopy}", "{snakemake.input.features}", "{args[pipeline_id]}"),
    file="{snakemake.output}", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE
)
__EOF

"""
)
