# -*- coding: utf-8 -*-
"""Render the calls of one VCF as a table for the report, with datavzrd.

The outputs are the calls table as TSV (one row per call) and the directory of the site, which the Snakemake report shows.
"""

from typing import TYPE_CHECKING

from snappy_wrappers.snappy_wrapper import PythonWrapper
from snappy_wrappers.tools.calls_report import Options, summarize
from snappy_wrappers.tools.qc_summary import render

if TYPE_CHECKING:
    from snakemake.iocontainers import snakemake

args = getattr(snakemake.params, "args", {})


def main():
    options = Options(
        page_size=args["page_size"],
        max_in_memory_rows=args["max_in_memory_rows"],
        linkouts=tuple(args["linkouts"]),
        assembly=args["assembly"],
    )
    summary = summarize(snakemake.input.vcf, options)
    render(summary, args["title"], snakemake.output.tsv, snakemake.output.html)


PythonWrapper(snakemake).run(main)
