# -*- coding: utf-8 -*-
"""The QC summary of a task: one job that renders the QC files of all its libraries as tables and
charts.

A step adds a subclass of :class:`QcSummaryStepPart` next to the step part that writes the QC
files. The job reads one set of files per library, writes the summary table as its result file
``qc_summary/out/qc_summary.tsv`` (one row per library) and renders a datavzrd site below
``qc_summary/report/html/``, which the Snakemake report shows (``snappy report``).
The ``qc_summary`` wrapper parses the files with the format module of
``snappy_wrappers.tools.qc_summary`` that ``format`` names.
"""

from __future__ import annotations

import itertools
import re
from functools import cached_property
from typing import Any

from snakemake.io import directory
from snakemake.iocontainers import Wildcards

from snappy_pipeline.utils import dictify
from snappy_pipeline.workflows.abstract import BaseStepPart, ResourceUsage


class QcSummaryStepPart(BaseStepPart):
    """Summary of the QC files of a task's libraries, rendered with datavzrd.

    A subclass names the ``source``: the step part and action whose outputs hold the QC files of
    one library; ``source_keys`` map a kind of file to the output key of that action, and
    ``format`` names the parser. ``report_outputs`` gives the subcategory, which titles the site,
    and the caption. ``entities`` lists the libraries that have the files, and the label of a
    file is the value of ``label_wildcard``, by default the entity.
    ``summary_options`` adds params that the format module takes as keyword arguments.
    """

    name = "qc_summary"
    actions = ("run",)

    #: Step part and action whose per-library outputs the summary reads
    source: tuple[str, str]
    #: Kind of file -> output key of the source action; the wrapper gets the files by kind
    source_keys: dict[str, str]
    #: Parser of the files: a module of ``snappy_wrappers.tools.qc_summary``
    format: str
    #: Wildcard of the source outputs whose value labels a file, else the entity
    label_wildcard: str | None = None

    result_keys = ("tsv", "tsv_md5")

    default_resource_usage = ResourceUsage(threads=1, runtime="2h", mem="4GB")

    def entities(self) -> list[str]:
        """The output entities that have the QC files: the task's, by default."""
        return self.parent.output_entities

    @cached_property
    def source_files(self) -> list[tuple[str, dict[str, str]]]:
        """``(label, path by kind)`` for each file set of each entity, from the output templates
        of the source action and the result wildcards of its entities."""
        part_name, action = self.source
        part = self.parent.sub_steps[part_name]
        outputs = part.get_output_files(action)
        templates = {
            kind: re.sub(r"\{(\w+),[^}]*\}", r"{\1}", str(outputs[key]))
            for kind, key in self.source_keys.items()
        }
        result = []
        for entity in self.entities():
            wildcards = part._result_wildcards(entity)
            names = list(wildcards)
            for values in itertools.product(*(wildcards[name] for name in names)):
                combination = dict(zip(names, values))
                label = combination[self.label_wildcard] if self.label_wildcard else entity
                files = {kind: tpl.format(**combination) for kind, tpl in templates.items()}
                result.append((label, files))
        return result

    @dictify
    def _get_input_files_run(self, wildcards: Wildcards):
        for kind in self.source_keys:
            yield kind, [files[kind] for _, files in self.source_files]

    @dictify
    def _get_output_files_run(self):
        tsv = "work/qc_summary/out/qc_summary.tsv"
        yield "tsv", tsv
        yield "tsv_md5", tsv + ".md5"
        yield "html", directory("work/qc_summary/report/html")

    def _get_log_file_run(self):
        return {"log": "work/qc_summary/log/qc_summary.log"}

    def summary_options(self) -> dict[str, Any]:
        """Further params for the format module's ``summarize``,
        such as the sample sheet facts a check compares the files with; none by default."""
        return {}

    def _get_params_run(self, wildcards: Wildcards) -> dict[str, Any]:
        return {
            "format": self.format,
            "labels": [label for label, _ in self.source_files],
            "title": self.report_outputs["run"]["html"].subcategory,
            **self.summary_options(),
        }

    def get_result_files(self) -> list[str]:
        """The summary table, if the task has libraries with QC files."""
        if not self.source_files:
            return []
        return ["output/qc_summary/out/qc_summary.tsv", "output/qc_summary/out/qc_summary.tsv.md5"]
