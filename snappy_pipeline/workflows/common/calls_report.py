# -*- coding: utf-8 -*-
"""The calls report of a task: one job per entity that renders the calls of its VCF as a table.

A step that writes VCFs registers :class:`CallsReportStepPart` when its ``calls_report`` section is set.
The job reads the task's own output VCF of the entity, writes the calls table as the result file ``{library_name}/report/calls/{library_name}.calls.tsv`` and renders a datavzrd site below ``{library_name}/report/calls/html/``, which the Snakemake report shows (``snappy report``).
The ``calls_report`` wrapper builds the table with ``snappy_wrappers.tools.calls_report``; the columns follow what the VCF carries.
"""

from __future__ import annotations

from typing import Any

from snakemake.io import directory
from snakemake.iocontainers import Wildcards

from snappy_pipeline.utils import dictify
from snappy_pipeline.workflows.abstract import BaseStepPart, ReportOutput, ResourceUsage


class CallsReportStepPart(BaseStepPart):
    """The calls of each VCF of the task as a table and a datavzrd site for the report."""

    name = "calls_report"
    actions = ("run",)
    entity_wildcard = "library_name"
    result_keys = ("tsv", "tsv_md5")

    default_resource_usage = ResourceUsage(threads=1, runtime="2h", mem="8GB")

    #: The site, and the table itself for download
    report_outputs = {
        "run": {
            "html": ReportOutput("calls", caption="calls_report.rst", htmlindex="index.html"),
            "tsv": ReportOutput("calls", caption="calls_report.rst"),
        }
    }

    def _get_input_files_run(self, wildcards: Wildcards) -> dict[str, str]:
        """The task's own output VCF of the entity (the PASS calls of mutect2)."""
        paths = self.parent.get_output_paths(
            self.parent.config, library_name=wildcards.library_name
        )
        return {"vcf": paths["vcf"], "vcf_tbi": paths["vcf_tbi"]}

    @dictify
    def _get_output_files_run(self):
        tsv = "work/{library_name}/report/calls/{library_name}.calls.tsv"
        yield "tsv", tsv
        yield "tsv_md5", tsv + ".md5"
        yield "html", directory("work/{library_name}/report/calls/html")

    def _get_log_file_run(self):
        return {"log": "work/{library_name}/log/{library_name}.calls_report.log"}

    def _get_params_run(self, wildcards: Wildcards) -> dict[str, Any]:
        settings = self.parent.config.calls_report
        return {
            "title": f"{self.parent.task_name}: calls of {wildcards.library_name}",
            "page_size": settings.page_size,
            "max_in_memory_rows": settings.max_in_memory_rows,
        }
