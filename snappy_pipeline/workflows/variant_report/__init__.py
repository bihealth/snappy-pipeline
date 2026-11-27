# -*- coding: utf-8 -*-
"""Report of the variant calls of a task

The step renders the calls of each VCF of ``depends_on.variants`` as a table.
The ``tool`` field is ``datavzrd``, which renders the table as a site for the report (``snappy report``).
The task writes one report per entity of the upstream task: per library, or per pedigree or donor when that task writes one VCF per related group.
The columns follow what the VCF carries:
varlociraptor's event probabilities and observations, or the genotypes, allele frequencies and depths of other callers;
the gene, consequence, impact and HGVS of annotated variants;
the gnomAD frequency and ClinVar classification that mehari adds.
The linkouts that depend on the genome build need an annotation task in the chain, which tags its variants with the build.
With ``igv`` set, the job also reads the BAMs of the entity's libraries (``depends_on.alignments``) and the genome (``depends_on.reference``),
and renders IGV views of the top calls of each view into the site.

Outputs
-------

Below ``tasks/<task>/output/``:

- ``{library_name}/report/calls/{library_name}.calls.tsv`` with ``.md5``, the calls of the VCF as a table (one row per call).
  The same table is a datavzrd site below ``report/calls/html/``, which ``snappy report`` shows.

The ``calls_report`` wrapper builds the table with ``snappy_wrappers.tools.calls_report``.
"""

from __future__ import annotations

from typing import Any

from snakemake.io import directory
from snakemake.iocontainers import Wildcards

from snappy_pipeline.utils import dictify
from snappy_pipeline.workflows.abstract import BaseStep, BaseStepPart, ReportOutput, ResourceUsage
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType, ExpectedVariants
from snappy_pipeline.workflows.variant_calling.scenario import event_suffix

from .model import Tool
from .model import VariantReport as VariantReportConfigModel

#: The genome builds that an annotation task tags its variants with
ASSEMBLIES = frozenset({"GRCh37", "GRCh38"})


class DatavzrdStepPart(BaseStepPart):
    """The calls of each VCF of the upstream task as a table and a datavzrd site for the report."""

    name = "datavzrd"
    actions = ("run",)
    wrappers = {"run": "calls_report"}
    result_keys = ("tsv", "tsv_md5")

    #: The site, and the table itself for download
    report_outputs = {
        "run": {
            "html": ReportOutput("calls", caption="calls_report.rst", htmlindex="index.html"),
            "tsv": ReportOutput("calls", caption="calls_report.rst"),
        }
    }

    default_resource_usage = ResourceUsage(threads=1, runtime="2h", mem="8GB")

    def tool_section(self):
        """datavzrd has no section of its own: the ``extra_arguments`` are the step's"""
        return self.config

    def _get_input_files_run(self, wildcards: Wildcards) -> dict[str, Any]:
        """The VCF of the entity from ``depends_on.variants``;
        with IGV views the BAMs of its libraries and the genome."""
        variants: ExpectedVariants = self.parent.get_upstream_paths(
            "variants", library_name=wildcards.library_name
        )
        inputs: dict[str, Any] = {"vcf": variants.vcf, "vcf_tbi": variants.vcf_tbi}
        if self.config.igv:
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

    def _get_params_run(self, wildcards: Wildcards) -> dict[str, Any]:
        return {
            "title": f"{self.parent.task_name}: calls of {wildcards.library_name}",
            "page_size": self.config.page_size,
            "max_in_memory_rows": self.config.max_in_memory_rows,
            "linkouts": [str(link) for link in self.config.linkouts],
            "columns": self.config.columns,
            "assembly": self.assembly(),
            "igv": self._igv_params(wildcards.library_name),
            "pedigree": self.pedigree(wildcards.library_name),
        }

    def igv_libraries(self, entity: str) -> list[str]:
        """The libraries whose BAMs the IGV views show:
        the members of the entity and, on cancer sheets, their matched normals."""
        libraries = list(self.parent.member_libraries(entity))
        for library in list(libraries):
            normal = self.parent.matched_normal(library)
            if normal and normal not in libraries:
                libraries.append(normal)
        return libraries

    def assembly(self) -> str | None:
        """The genome build the variants are tagged with (``variant_annotation`` adds it and the steps after it keep it), else ``None``."""
        tags = self.parent.upstream_signature("variants").tags
        return next((tag for tag in tags if tag in ASSEMBLIES), None)

    def pedigree(self, entity: str) -> list[dict[str, Any]]:
        """The members of the entity's pedigree for the pedigree plot: every donor of the
        sample sheet's pedigree, named after its library (the VCF sample) or, unsequenced,
        after the donor, with its parents, sex and affected status, and for a sequenced one
        its part of a per-member event name (``denovo_<member>``), so the plot marks the member
        of such an event; empty for an entity that is no pedigree."""
        if self.parent.effective_group_by == "library":
            return []
        donors = self.parent.pedigree_donors(entity, include_background=True)
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
                    "event_suffix": event_suffix(row._asdict()) if row.library_name else None,
                }
            )
        return members

    def _igv_params(self, entity: str) -> dict[str, Any] | None:
        igv = self.config.igv
        if not igv:
            return None
        return {
            "top": igv.top,
            "views": list(igv.views),
            "flanking": igv.flanking,
            "max_span": igv.max_span,
            "max_size_mb": igv.max_size_mb,
            "extra_arguments": list(igv.extra_arguments),
            "libraries": self.igv_libraries(entity),
        }


class VariantReportWorkflow(BaseStep):
    """Render the calls of the VCFs of a task as tables for the report"""

    name = "variant_report"

    produces = [DataSignature(DataType.TABULAR, frozenset({"calls"}))]
    config_model_class = VariantReportConfigModel
    tool_parts = {Tool.datavzrd: DatavzrdStepPart}

    @classmethod
    def get_output_paths(cls, config, signature=None, **kwargs) -> dict[str, str]:
        """The calls table of an entity, a ``work/`` path."""
        lib = kwargs.get("library_name", "{library_name}")
        return {"tsv": f"work/{lib}/report/calls/{lib}.calls.tsv"}
