# -*- coding: utf-8 -*-
"""Scatter-gather over regions of the genome.

A step part that processes the genome in parallel inherits :class:`RegionScatter` before ``BaseStepPart``.
Its ``scatter`` action splits the reference's ``.fai`` into the ``scattergather`` parts named after the step part, with the ``split_regions`` wrapper.
The regions are written once per task below ``work/scatter/``, for all entities of the task.
The step Snakefile declares ``scattergather: <part name>=<number of parts>``.
"""

from __future__ import annotations

from typing import Any, Callable

from snappy_pipeline.workflows.abstract import ResourceUsage


class ScatterGather:
    """Snakemake's ``scatter`` and ``gather`` functions of a step part, named after it."""

    def _scattergather(self, kind: str) -> Callable[..., Any]:
        """``scatter.<name>`` or ``gather.<name>`` of the Snakefile's ``scattergather``"""
        return getattr(self.parent.workflow.globals.get(kind), self.name)

    def gathered(self, template: str, wildcards) -> list[str]:
        """All parts of ``template``: ``{scatteritem}`` expanded by ``gather``, the other wildcards filled in."""
        partial = template.replace("{scatteritem}", "{{scatteritem}}").format(**wildcards)
        return self._scattergather("gather")(partial)


class RegionScatter(ScatterGather):
    """The ``scatter`` action of a step part: the regions of the genome, one BED file per part."""

    #: The regions of one part, shared by the entities of the task
    region_tpl = "work/scatter/{scatteritem}.region.bed"

    #: The wrapper of the ``scatter`` action, for the step part's ``wrappers``
    scatter_wrapper = {"scatter": "split_regions"}

    #: Bases added on both sides of each region; 0 keeps the regions disjoint
    scatter_padding = 0

    scatter_resource_usage = ResourceUsage(threads=1, runtime="2m", mem="1000MB")

    def scatter_ignore_contigs(self) -> list[str]:
        """The regular expressions of the contigs that the regions leave out"""
        raise NotImplementedError

    def get_resource_usage(self, action: str) -> ResourceUsage:
        if action == "scatter":
            return self.scatter_resource_usage
        return super().get_resource_usage(action)

    def _get_input_files_scatter(self, wildcards) -> dict[str, Any]:
        return {"fai": self.parent.get_upstream_paths("reference").fai}

    def _get_params_scatter(self, wildcards) -> dict[str, Any]:
        return {"ignore_contigs": self.scatter_ignore_contigs(), "padding": self.scatter_padding}

    def _get_output_files_scatter(self) -> dict[str, Any]:
        return {"regions": self._scattergather("scatter")(self.region_tpl)}

    def _get_log_file_scatter(self) -> dict[str, str]:
        return {"log": "work/scatter/log.log"}

    def region(self, wildcards) -> str:
        """The regions of the part of a job (its ``scatteritem``)"""
        return self.region_tpl.format(**wildcards)
