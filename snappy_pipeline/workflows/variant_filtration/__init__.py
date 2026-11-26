# -*- coding: utf-8 -*-
"""Variant filtration

The step filters or tags the variants of the task named in ``depends_on.variants``.
The ``tool`` field selects the filter, and one task applies one filter.
Chain several tasks for a sequence of filters.

``bcftools``
    Filtering by a ``bcftools`` expression in ``include`` or ``exclude``.

``vembrane``
    Expression-based filtering with ``vembrane``: Python expressions on the fields of a record, with ``ANN`` parsed by default.
    It is the recommended filter for genotype-based criteria, such as inheritance patterns in a family.
    With ``group_by: related`` upstream, one VCF holds the whole pedigree, and an expression addresses samples by name.
    For example, ``is_het("child-N1-DNA1-WGS1") and is_hom_ref("father-N1-DNA1-WGS1") and is_hom_ref("mother-N1-DNA1-WGS1")`` selects de novo variants.
    Helpers that several expressions share go into ``context_files``.

``regions``
    Filtering by the regions of a BED file, from the task in ``depends_on.regions``.
    ``action: exclude`` applies the filter to the variants inside the regions, ``include`` to those outside.

``dkfz``
    DKFZ bias filter.
    It needs ``depends_on.alignments`` and ``depends_on.reference``.

``ebfilter``
    EBFilter.
    Variants with one of the ``skip_filters`` skip EBFilter, if the VCF header declares the filter.
    By default, ``skip_filters`` holds the ``germline``, ``weak_evidence`` and ``OffExome`` filters of Mutect2.
    With an annotated input, the variants whose annotation contains ``skip_annotation`` skip EBFilter too.
    By default, ``skip_annotation`` is ``stream_gene_variant``, which marks the variants in or near genes.
    It needs ``depends_on.alignments`` and ``depends_on.reference``.

``varlociraptor``
    FDR control of the calls of a ``variant_calling`` task with ``tool: varlociraptor``.
    ``varlociraptor filter-calls control-fdr`` selects the records whose posterior probabilities for the ``events`` meet the ``fdr`` in the chosen ``fdr_mode``.
    The control runs once per variant type.
    Chain this task right after the calling task, before annotation.
    The events are those of the calling task's scenario, such as ``somatic_tumor_high`` and ``somatic_tumor_medium``.
    The germline events are ``present``, ``denovo_*`` (a pattern for the per-member events ``denovo_<member>``), ``recessive`` and ``dominant``.
    The germline events partition the carried variants, so a record passes only for the events listed.
    In smart mode, a record passes when the listed events hold at least half of its probability.
    List all four germline events to keep every variant.
    An entry of ``events`` is an event name, which must exist in the calls, or a shell-style pattern, which may match none.
    The names and patterns work the same for the events of a scenario ``template``.
    The ``mode`` of ``varlociraptor`` defaults to ``filter``, so the task after it reads a call set.

Every filter has a ``mode``: ``tag`` marks variants in ``FILTER``, ``filter`` removes them.
``bcftools``, ``regions``, ``ebfilter`` and ``varlociraptor`` in ``tag`` mode use the task name as ``FILTER`` ID, so each task of a chain has its own.
``vembrane`` uses the names of its ``expressions``, and ``dkfz`` the IDs of the DKFZ bias filter (``bPcr``, ``bSeq``).
No other task of the project may be named like one of these IDs or like ``PASS``.
No task may write an ID that a task upstream of it writes, such as an expression named after a filter of mutect2.
The task produces one VCF per VCF of its upstream task (per tumor library, or per pedigree for germline calls).
Its output has the input's tags plus ``filtered``.

The example tags the variants in regions of low mappability, from a BED file of an ``external_data`` task:

.. code-block:: yaml

    tasks:
      - step: external_data
        name: calls
        config:
          produces: {type: variants, tags: [somatic, snv, indel]}
          library_files: [calls.tsv]   # the columns library, vcf, vcf_tbi

      - step: external_data
        name: low_mappability
        config:
          produces: {type: regions, tags: [mappability]}
          files: {bed: /refs/GRCh38/low_mappability.bed}

      - step: variant_filtration
        name: mappability_filter
        config:
          depends_on:
            variants: calls
            regions: low_mappability
          tool: regions
          regions: {action: exclude, mode: tag}

Outputs
-------

Below ``tasks/<task>/output/``:

- ``{library_name}/out/{library_name}.vcf.gz`` with ``.tbi``, ``.md5`` and ``.tbi.md5`` files.

Downstream tasks read ``vcf`` and ``vcf_tbi`` through ``depends_on``.
"""

import random
from functools import cached_property
from typing import Any

from snakemake.iocontainers import Wildcards

from snappy_pipeline.base import ConfigError
from snappy_pipeline.models import ToolModel
from snappy_pipeline.utils import dictify
from snappy_pipeline.workflows.abstract import (
    MATCHED_NORMAL_RELATIONSHIP,
    BaseStep,
    BaseStepPart,
    ResourceUsage,
)
from snappy_pipeline.workflows.common import vcf_outputs
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType

from .model import Ebfilter as EbfilterConfig
from .model import Tool
from .model import VariantFiltration as VariantFiltrationConfigModel

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

#: The filtered VCF of an entity, without ``.vcf.gz``
_OUT_PREFIX = "work/{library_name}/out/{library_name}"


# ---------------------------------------------------------------------------
# Base step part
# ---------------------------------------------------------------------------
class VariantFiltrationStepPart(BaseStepPart):
    """Base class for all variant filtration tool step parts.

    All concrete tool step parts use ``name = "filter"`` because only one is
    ever registered per workflow task instance.
    """

    name = "filter"
    actions = ("run",)
    entity_wildcard = "library_name"
    result_keys = {"run": tuple(vcf_outputs(_OUT_PREFIX))}

    resource_usage = {
        "run": ResourceUsage(threads=1, runtime="4h", mem=f"{8 * 1024}MB"),
    }

    def tool_section(self) -> ToolModel:
        """The config section of the task's tool; all filter step parts are named ``filter``"""
        return getattr(self.config, self.config.tool)

    def _get_log_file(self, action):
        return {"log": "work/{library_name}/log/{library_name}.log"}

    @dictify
    def _get_input_files_run(self, wildcards: Wildcards):
        variant = self.parent.get_upstream_paths("variants", library_name=wildcards.library_name)
        yield "vcf", variant.vcf
        yield "vcf_tbi", variant.vcf_tbi

    def _get_output_files_run(self):
        return vcf_outputs(_OUT_PREFIX)

    def _get_params_run(self, wildcards: Wildcards) -> dict[str, Any]:
        # The task name keeps FILTER IDs unique along a chain of filter tasks.
        return {"filter_name": self.parent.task_name} | self.tool_section().model_dump(
            by_alias=True
        )


# ---------------------------------------------------------------------------
# Concrete tool step parts
# ---------------------------------------------------------------------------
class BcftoolsStepPart(VariantFiltrationStepPart):
    """bcftools expression filter."""

    wrappers = {"run": "bcftools/filter"}


class VembraneStepPart(VariantFiltrationStepPart):
    """Unified vembrane tag/filter step."""

    wrappers = {"run": "vembrane/tag"}

    @dictify
    def _get_input_files_run(self, wildcards: Wildcards):
        yield from super()._get_input_files_run(wildcards).items()
        cfg = self.config.vembrane
        if cfg.aux:
            yield "aux", list(cfg.aux.values())
        if cfg.context_files:
            yield "context_files", list(cfg.context_files)
        if cfg.ontology:
            yield "ontology", cfg.ontology

    def _get_params_run(self, wildcards: Wildcards) -> dict[str, Any]:
        cfg = self.config.vembrane
        # vembrane tag and filter both take the aux files, the context and the ontology
        params: dict[str, Any] = {
            "mode": cfg.mode,
            "aux_names": list(cfg.aux),
            "context": cfg.context,
        }
        if cfg.mode == "tag":
            params["expressions"] = cfg.expressions
        else:
            params["expression"] = cfg.expression
        return params


class RegionsStepPart(VariantFiltrationStepPart):
    """Region/BED-based filter via bcftools."""

    wrappers = {"run": "bcftools/regions"}

    @dictify
    def _get_input_files_run(self, wildcards: Wildcards):
        yield from super()._get_input_files_run(wildcards).items()
        # The wrapper reads the BED file as input "include" or "exclude"
        yield self.config.regions.action, self.parent.get_upstream_paths("regions").bed

    def _get_params_run(self, wildcards: Wildcards) -> dict[str, Any]:
        return {"filter_name": self.parent.task_name, "mode": self.config.regions.mode}


# ---------------------------------------------------------------------------
# BAM-aware base (dkfz / ebfilter)
# ---------------------------------------------------------------------------
class _BamAwareStepPart(VariantFiltrationStepPart):
    """Mixin that adds tumor (and optionally normal) BAM inputs."""

    @dictify
    def _get_input_files_run(self, wildcards: Wildcards):
        yield from super()._get_input_files_run(wildcards).items()

        yield from self.parent.reference_inputs().items()

        lib = self.parent.primary_library(wildcards.library_name)
        yield "bam", self.parent.get_upstream_paths("alignments", library_name=lib).bam

        if normal_lib := self.parent.matched_normal(lib):
            normal = self.parent.get_upstream_paths("alignments", library_name=normal_lib)
            yield "normal", normal.bam


class DkfzStepPart(_BamAwareStepPart):
    """DKFZ bias filter."""

    wrappers = {"run": "dkfz_bias_filter"}

    resource_usage = {
        "run": ResourceUsage(threads=1, runtime="12h", mem=f"{3 * 1024}MB"),
    }

    def _get_params_run(self, wildcards: Wildcards) -> dict[str, Any]:
        # The filter writes its own FILTER IDs (bPcr, bSeq), not the task name
        params = super()._get_params_run(wildcards)
        params.pop("filter_name")
        return params


class EbfilterStepPart(_BamAwareStepPart):
    """EBFilter."""

    actions = ("run", "write_panel")
    wrappers = {"run": "eb_filter"}
    extra_arguments_action = "run"

    resource_usage = {
        "run": ResourceUsage(threads=1, runtime="24h", mem=f"{2 * 1024}MB"),
        "write_panel": ResourceUsage(threads=1, runtime="1h", mem=f"{2 * 1024}MB"),
    }

    @dictify
    def _get_input_files_run(self, wildcards: Wildcards):
        yield from super()._get_input_files_run(wildcards).items()
        cfg: EbfilterConfig = self.config.ebfilter
        panel_file = cfg.path_panel_of_normals_sample_list
        if not panel_file:
            panel_file = self._get_output_files_write_panel()["txt"]
        yield "txt", panel_file

    def _get_params_run(self, wildcards: Wildcards) -> dict[str, Any]:
        params = super()._get_params_run(wildcards)
        # The panel file and the options that choose it are inputs of the rule, not of EBFilter
        for key in ("path_panel_of_normals_sample_list", "panel_of_normals_size", "shuffle_seed"):
            params.pop(key, None)
        return params | {
            "has_annotation": self.parent.variants_are_annotated,
        }

    def _get_input_files_write_panel(self, wildcards):
        return {}

    def _get_output_files_write_panel(self):
        return {"txt": "work/panel_of_normals/out/panel_of_normals.txt"}

    def _get_params_write_panel(self, wildcards):
        """The panel options and the libraries they choose, so that a change writes the panel again"""
        cfg: EbfilterConfig = self.config.ebfilter
        return {
            "panel_of_normals_size": cfg.panel_of_normals_size,
            "shuffle_seed": cfg.shuffle_seed,
            "libraries": self._panel_of_normal_libraries(),
        }

    def write_panel_of_normals_file(self, wildcards):
        out_path = self._get_output_files_write_panel()["txt"]
        with open(out_path, "wt") as fh:
            for bam in self._get_panel_of_normal_bams(wildcards):
                print(bam, file=fh)

    def _panel_of_normal_libraries(self) -> list[str]:
        """The normal DNA libraries of the panel: a seeded shuffle, cut to the panel size"""
        df = self.parent.build_library_dataframe()
        normal_df = df[(df["role"] == "normal") & (df["extraction_type"] == "dna")]
        libraries = sorted(normal_df["library_name"].tolist())

        cfg: EbfilterConfig = self.config.ebfilter
        random.seed(cfg.shuffle_seed)
        random.shuffle(libraries)
        return libraries[: cfg.panel_of_normals_size]

    def _get_panel_of_normal_bams(self, wildcards) -> list[str]:
        return [
            self.parent.get_upstream_paths("alignments", library_name=name).bam
            for name in self._panel_of_normal_libraries()
        ]


class VarlociraptorStepPart(VariantFiltrationStepPart):
    """FDR control of varlociraptor calls."""

    wrappers = {"run": "varlociraptor/control_fdr"}

    def _get_params_run(self, wildcards: Wildcards) -> dict[str, Any]:
        # Only tag mode writes a FILTER ID, the task name
        params = super()._get_params_run(wildcards)
        if params["mode"] == "filter":
            params.pop("filter_name")
        return params


# ---------------------------------------------------------------------------
# Workflow
# ---------------------------------------------------------------------------
class VariantFiltrationWorkflow(BaseStep):
    """Unified variant filtration - applies a single configurable filter tool to any VCF input.

    Supports somatic-specific tools (dkfz, ebfilter) as well as generic tools
    (bcftools, vembrane, regions). Use multiple tasks with
    ``depends_on.variants`` to compose a sequential filter pipeline.
    """

    name = "variant_filtration"
    produces = [DataSignature(DataType.VARIANTS, frozenset({"filtered"}))]
    config_model_class = VariantFiltrationConfigModel

    tool_parts = {
        Tool.bcftools: BcftoolsStepPart,
        Tool.vembrane: VembraneStepPart,
        Tool.regions: RegionsStepPart,
        Tool.dkfz: DkfzStepPart,
        Tool.ebfilter: EbfilterStepPart,
        Tool.varlociraptor: VarlociraptorStepPart,
    }

    default_relationships = MATCHED_NORMAL_RELATIONSHIP

    @classmethod
    def task_produces(cls, config, upstream):
        """The variants this task reads, tagged ``filtered``; the FDR control of varlociraptor
        needs varlociraptor calls (tag ``varlociraptor``, which annotation keeps)."""
        variants = upstream["variants"]
        if config.tool == Tool.varlociraptor and "varlociraptor" not in variants.tags:
            raise ConfigError(
                "tool: varlociraptor controls the FDR of the calls of a variant_calling task with "
                "tool: varlociraptor, but depends_on.variants produces variants tagged "
                + ", ".join(sorted(variants.tags))
            )
        return (variants.with_tags("filtered"),)

    @classmethod
    def get_output_paths(cls, config, signature=None, **kwargs) -> dict[str, str]:
        """Return local filtered-variant output paths for downstream consumers."""
        lib = kwargs.get("library_name", "{library_name}")
        paths = vcf_outputs(_OUT_PREFIX.format(library_name=lib))
        return {key: paths[key] for key in ("vcf", "vcf_tbi")}

    @classmethod
    def filter_ids(cls, config):
        """The IDs of the DKFZ bias filter, and the expression names of ``vembrane`` in tag mode."""
        if config.tool == Tool.dkfz:
            return frozenset({"bPcr", "bSeq"})
        if config.tool == Tool.vembrane and config.vembrane.mode == "tag":
            return frozenset(config.vembrane.expressions)
        return frozenset()

    @cached_property
    def variants_are_annotated(self) -> bool:
        """Whether the variants of ``depends_on.variants`` carry an annotation (tag ``annotated``)."""
        return "annotated" in self.upstream_signature("variants").tags
