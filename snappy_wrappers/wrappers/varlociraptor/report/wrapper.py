# -*- coding: utf-8 -*-
"""Render the varlociraptor calls of one VCF as a table for the report, with datavzrd.

The outputs are the calls table as TSV (one row per call, the event probabilities and the
per-sample values) and the directory of the site, which the Snakemake report shows.
"""

from typing import TYPE_CHECKING

from snappy_wrappers.snappy_wrapper import PythonWrapper
from snappy_wrappers.tools.qc_summary import render
from snappy_wrappers.tools.qc_summary.varlociraptor_calls import summarize

if TYPE_CHECKING:
    from snakemake.iocontainers import snakemake

args = getattr(snakemake.params, "args", {})


def main():
    summary = summarize(snakemake.input.vcf)
    render(summary, args["title"], snakemake.output.tsv, snakemake.output.html)


PythonWrapper(snakemake).run(main)
