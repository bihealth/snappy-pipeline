# -*- coding: utf-8 -*-
"""Render the calls of one VCF as a table for the report, with datavzrd.

The outputs are the calls table as TSV (one row per call) and the directory of the site,
which the Snakemake report shows.
With ``igv`` set,
igv-reports renders the top calls of the views into ``igv/top.html`` of the site first,
so that the table can link to it.
"""

import tempfile
from typing import TYPE_CHECKING

from snappy_wrappers.snappy_wrapper import PythonWrapper
from snappy_wrappers.tools.calls_report import Options, igv, summarize
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
        pedigree=tuple(args["pedigree"]),
    )
    summary = summarize(snakemake.input.vcf, options)
    settings = args["igv"]
    with tempfile.TemporaryDirectory() as workdir:
        page = None
        if settings:
            page, sites, requested = igv.build_page(
                summary,
                settings,
                snakemake.input.bam,
                snakemake.input.bai,
                snakemake.input.vcf,
                snakemake.input.vcf_tbi,
                snakemake.input.reference,
                workdir,
                args["title"],
            )
            igv.annotate(summary, sites, requested, settings["max_size_mb"])
        render(summary, args["title"], snakemake.output.tsv, snakemake.output.html)
        igv.install_page(page, snakemake.output.html)


PythonWrapper(snakemake).run(main)
