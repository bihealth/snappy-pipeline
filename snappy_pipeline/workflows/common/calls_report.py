# -*- coding: utf-8 -*-
"""The calls report of a task: one job per entity that renders the calls of its VCF as a table.

A step that writes VCFs registers :class:`CallsReportStepPart` when its ``calls_report`` section is set.
With ``calls_report.igv`` the job also reads the BAMs of the entity's libraries (``depends_on.alignments``)
and the genome, and writes IGV views of the top calls into the site.
The job reads the task's own output VCF of the entity,
writes the calls table as the result file ``{library_name}/report/calls/{library_name}.calls.tsv`` and renders a datavzrd site below ``{library_name}/report/calls/html/``,
which the Snakemake report shows (``snappy report``).
The ``calls_report`` wrapper builds the table with ``snappy_wrappers.tools.calls_report``;
the columns follow what the VCF carries.
"""

from __future__ import annotations

from typing import Any

from snakemake.io import directory
from snakemake.iocontainers import Wildcards

from snappy_pipeline.utils import dictify
from snappy_pipeline.workflows.abstract import BaseStepPart, ReportOutput, ResourceUsage

#: The genome builds that an annotation task tags its variants with
ASSEMBLIES = frozenset({"GRCh37", "GRCh38"})


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

    def igv_libraries(self, entity: str) -> list[str]:
        """The libraries whose BAMs the IGV views show:
        the members of the entity and, on cancer sheets, their matched normals."""
        libraries = list(self.parent.member_libraries(entity))
        table = self.parent.libraries
        if "matched_normal_lib" in table.columns:
            for library in list(libraries):
                rows = table[table["library_name"] == library]
                normal = rows.iloc[0].get("matched_normal_lib") if not rows.empty else None
                if normal and normal not in libraries:
                    libraries.append(normal)
        return libraries

    def _get_input_files_run(self, wildcards: Wildcards) -> dict[str, Any]:
        """The task's own output VCF of the entity (the PASS calls of mutect2);
        with IGV views the BAMs of its libraries and the genome."""
        paths = self.parent.get_output_paths(
            self.parent.config, library_name=wildcards.library_name
        )
        inputs: dict[str, Any] = {"vcf": paths["vcf"], "vcf_tbi": paths["vcf_tbi"]}
        if self.parent.config.calls_report.igv:
            alignments = [
                self.parent.get_upstream_paths("alignments", library_name=library)
                for library in self.igv_libraries(wildcards.library_name)
            ]
            inputs["bam"] = [paths.bam for paths in alignments]
            inputs["bai"] = [paths.bai for paths in alignments]
            inputs |= self.parent.reference_inputs()
        return inputs

    @dictify
    def _get_output_files_run(self):
        tsv = "work/{library_name}/report/calls/{library_name}.calls.tsv"
        yield "tsv", tsv
        yield "tsv_md5", tsv + ".md5"
        yield "html", directory("work/{library_name}/report/calls/html")

    def _get_log_file_run(self):
        return {"log": "work/{library_name}/log/{library_name}.calls_report.log"}

    def assembly(self) -> str | None:
        """The genome build the task's variants are tagged with (``variant_annotation`` adds it and the steps after it keep it), else ``None``."""
        for signature in self.parent.project.signatures[self.parent.task_name]:
            for tag in signature.tags:
                if tag in ASSEMBLIES:
                    return tag
        return None

    def _get_params_run(self, wildcards: Wildcards) -> dict[str, Any]:
        settings = self.parent.config.calls_report
        return {
            "title": f"{self.parent.task_name}: calls of {wildcards.library_name}",
            "page_size": settings.page_size,
            "max_in_memory_rows": settings.max_in_memory_rows,
            "linkouts": [str(link) for link in settings.linkouts],
            "assembly": self.assembly(),
            "igv": self._igv_params(wildcards.library_name),
            "pedigree": self.pedigree(wildcards.library_name),
        }

    def pedigree(self, entity: str) -> list[dict[str, Any]]:
        """The members of the entity's pedigree for the pedigree plot: every donor of the
        sample sheet's pedigree, named after its library (the VCF sample) or, unsequenced,
        after the donor, with its parents, sex and affected status; empty for an entity that
        is no pedigree."""
        if self.parent.effective_group_by != "related":
            return []
        donors = self.parent.pedigree_donors(entity)
        if donors.empty:
            return []
        table = self.parent.build_library_dataframe()
        status = dict(zip(table["library_name"], table["disease_state"]))
        names = {
            row.donor_name: row.library_name or row.donor_name
            for row in donors.itertuples(index=False)
        }
        members = []
        for row in donors.itertuples(index=False):
            affected = {2: True, 1: False}.get(status.get(row.library_name))
            members.append(
                {
                    "name": names[row.donor_name],
                    "label": row.donor_name,
                    "father": names.get(row.father_name, "0"),
                    "mother": names.get(row.mother_name, "0"),
                    "sex": row.sex,
                    "affected": affected,
                    "sample": row.library_name or None,
                }
            )
        return members

    def _igv_params(self, entity: str) -> dict[str, Any] | None:
        igv = self.parent.config.calls_report.igv
        if not igv:
            return None
        return {
            "top": igv.top,
            "views": list(igv.views),
            "flanking": igv.flanking,
            "max_size_mb": igv.max_size_mb,
            "extra_arguments": list(igv.extra_arguments),
            "libraries": self.igv_libraries(entity),
        }
