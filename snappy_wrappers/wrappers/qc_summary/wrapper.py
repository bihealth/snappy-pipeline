# -*- coding: utf-8 -*-
"""Render the QC files of a task's libraries as tables and charts with datavzrd.

The inputs are the files by kind, one per library in the order of ``args["labels"]``;
``args["format"]`` names the parser (``snappy_wrappers.tools.qc_summary``),
and the other params are its options (the sample sheet facts of a check).
The outputs are the summary TSV (one row per library) and the directory of the site,
which the Snakemake report shows.
"""

from typing import TYPE_CHECKING

from snappy_wrappers.snappy_wrapper import PythonWrapper
from snappy_wrappers.tools.qc_summary import render, summarize

if TYPE_CHECKING:
    from snakemake.iocontainers import snakemake

args = getattr(snakemake.params, "args", {})


def main():
    files = {kind: list(paths) for kind, paths in snakemake.input.items()}
    options = {
        key: value for key, value in args.items() if key not in ("format", "labels", "title")
    }
    summary = summarize(args["format"], files, list(args["labels"]), **options)
    render(summary, args["title"], snakemake.output.tsv, snakemake.output.html)


PythonWrapper(snakemake).run(main)
