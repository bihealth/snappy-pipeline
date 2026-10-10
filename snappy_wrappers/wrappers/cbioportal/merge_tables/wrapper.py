# -*- coding: utf-8 -*-
"""Wrapper for merging multiple tables in R on shared columns"""

import os

from snappy_wrappers.snappy_wrapper import ShellWrapper

__author__ = "Eric Blanc <eric.blanc@bih-charite.de>"

r_script = os.path.abspath(os.path.join(os.path.dirname(__file__), "script.R"))
helper_functions = os.path.join(os.path.dirname(r_script), "..", "helper_functions.R")

args = getattr(snakemake.params, "args", {})

tables = zip(args["samples"], snakemake.input.tables, strict=True)
filenames = ", ".join(['"{}"="{}"'.format(sample, table) for sample, table in tables])
mappings = snakemake.input.get("mappings", "")
extra = dict(args.get("extra_args", {}))
if "features" in snakemake.input.keys():
    extra["tx_obj"] = snakemake.input.features
extra_args = ", ".join(['"{}"="{}"'.format(str(k), str(v)) for k, v in extra.items()])

ShellWrapper(snakemake).run(
    r"""
# Run the R script --------------------------------------------------------------------------------

R --vanilla --slave << __EOF
source("{helper_functions}")
source("{r_script}")
write.table(
    merge_tables(list({filenames}), mappings="{mappings}", type="{args[action_type]}", args=list({extra_args})),
    file="{snakemake.output}", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE
)
__EOF
"""
)
