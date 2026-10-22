# -*- coding: utf-8 -*-
"""Variant filtration

Filters or tags the variants of the task named in ``depends_on.variants``.
The ``tool`` field selects the filter, and one task applies one filter.
Chain several tasks for a sequence of filters.

``bcftools``
    Expression-based filtering with ``include`` or ``exclude``.

``vembrane``
    Expression-based filtering with ``vembrane``.

``regions``
    Filtering by the regions of a BED file, from the task in ``depends_on.regions``.
    ``action: exclude`` applies the filter to the variants inside the regions, ``include`` to those outside.

``dkfz``
    DKFZ bias filter.
    Needs ``depends_on.alignments`` and ``depends_on.reference``.

``ebfilter``
    EBFilter.
    Needs ``depends_on.alignments`` and ``depends_on.reference``.

Every filter has a ``mode``: ``tag`` marks variants in ``FILTER``, ``filter`` removes them.
``bcftools``, ``regions`` and ``ebfilter`` use the task name as ``FILTER`` ID, so each task of a chain has its own; ``vembrane`` uses the names of its ``expressions``, ``dkfz`` the IDs of the DKFZ bias filter (``bPcr``, ``bSeq``).
No other task of the project may be named like one of these or like ``PASS``, and no task may write an ID that a task upstream of it writes (such as an expression named after a filter of mutect2).
The task produces one VCF per tumor library.
Its output has the input's tags plus ``filtered``.

The example tags the variants in regions of low mappability, from a BED file of an ``external_data`` task:

.. code-block:: yaml

    tasks:
      - step: external_data
        name: calls
        config:
          produces: {type: variants, tags: [somatic, snv, indel]}
          search_paths: [/data/calls]
          search_patterns:
            - {vcf: '.+\\.vcf\\.gz', vcf_tbi: '.+\\.vcf\\.gz\\.tbi'}

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
from typing import Any

from snakemake.iocontainers import Wildcards

from snappy_pipeline.models import RelationshipDefinition
from snappy_pipeline.utils import dictify, listify
from snappy_pipeline.workflows.abstract import (
    BaseStep,
    BaseStepPart,
    ResourceUsage,
)
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType
from snappy_pipeline.workflows.ngs_mapping.model import ExpectedAlignments

from .model import (
    Ebfilter as EbfilterConfig,
)
from .model import (
    VariantFiltration as VariantFiltrationConfigModel,
)

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

#: Output keys and extensions of the filtered VCF
_OUTPUTS = {
    "vcf": ".vcf.gz",
    "vcf_tbi": ".vcf.gz.tbi",
    "vcf_md5": ".vcf.gz.md5",
    "vcf_tbi_md5": ".vcf.gz.tbi.md5",
}


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
    resource_usage = {
        "run": ResourceUsage(threads=1, runtime="4h", mem=f"{8 * 1024}MB"),
    }
    published_keys = {"run": tuple(_OUTPUTS)}
    entity_wildcard = "tumor_library"

    @dictify
    def _get_input_files_run(self, wildcards: Wildcards):
        lib = wildcards.tumor_library
        variant = self.parent.get_upstream_paths("variants", library_name=lib)
        # Accept both typed schema and plain dict
        if isinstance(variant, dict):
            yield "vcf", variant["vcf"]
            yield "vcf_tbi", variant["vcf_tbi"]
        else:
            yield "vcf", variant.vcf
            yield "vcf_tbi", variant.vcf_tbi

    def _get_output_files_run(self):
        return {
            key: f"work/{{tumor_library}}/out/{{tumor_library}}{ext}"
            for key, ext in _OUTPUTS.items()
        }

    def _get_log_file(self, action):
        return {"log": "work/{tumor_library}/log/{tumor_library}.log"}

    def _get_params_run(self, wildcards: Wildcards) -> dict[str, Any]:
        cfg = getattr(self.config, self.config.tool)
        params: dict[str, Any] = {
            # The task name keeps FILTER IDs unique along a chain of filter tasks.
            "filter_name": self.parent.task_name,
        }
        if cfg is not None:
            params.update(cfg.model_dump(by_alias=True))
        return params


# ---------------------------------------------------------------------------
# BAM-aware base (dkfz / ebfilter)
# ---------------------------------------------------------------------------


class _BamAwareStepPart(VariantFiltrationStepPart):
    """Mixin that adds tumor (and optionally normal) BAM inputs."""

    def __init__(self, parent):
        super().__init__(parent)

    @dictify
    def _get_input_files_run(self, wildcards: Wildcards):
        yield from super()._get_input_files_run(wildcards).items()

        yield from self.parent.reference_inputs().items()

        lib = wildcards.tumor_library
        tumor_aln: ExpectedAlignments = self.parent.get_upstream_paths(
            "alignments", library_name=lib
        )
        if isinstance(tumor_aln, dict):
            yield "bam", tumor_aln["bam"]
        else:
            yield "bam", tumor_aln.bam

        df = self.parent.build_library_dataframe()
        tumor_df = df[df["library_name"] == lib]
        normal_lib = None
        if not tumor_df.empty:
            normal_lib = tumor_df.iloc[0].get("matched_normal_lib") or None
        if normal_lib:
            normal_aln: ExpectedAlignments = self.parent.get_upstream_paths(
                "alignments", library_name=normal_lib
            )
            if isinstance(normal_aln, dict):
                yield "normal", normal_aln["bam"]
            else:
                yield "normal", normal_aln.bam


# ---------------------------------------------------------------------------
# Concrete tool step parts
# ---------------------------------------------------------------------------


class BcftoolsStepPart(VariantFiltrationStepPart):
    """bcftools expression filter."""


class VembraneStepPart(VariantFiltrationStepPart):
    """Unified vembrane tag/filter step."""

    def _get_params_run(self, wildcards: Wildcards) -> dict[str, Any]:
        cfg = self.config.vembrane
        if cfg is None:
            raise ValueError("vembrane configuration is required")
        params: dict[str, Any] = {"mode": cfg.mode, "extra_args": cfg.extra_args}
        if cfg.mode == "tag":
            params["expressions"] = cfg.expressions
        else:
            params.update(
                {
                    "expression": cfg.expression,
                    "aux_names": list(cfg.aux),
                    "context": cfg.context,
                }
            )
        return params

    @dictify
    def _get_input_files_run(self, wildcards: Wildcards):
        yield from super()._get_input_files_run(wildcards).items()
        cfg = self.config.vembrane
        if cfg is None or cfg.mode != "filter":
            return
        if cfg.aux:
            yield "aux", list(cfg.aux.values())
        if cfg.context_files:
            yield "context_files", list(cfg.context_files)
        if cfg.ontology:
            yield "ontology", cfg.ontology


class RegionsStepPart(VariantFiltrationStepPart):
    """Region/BED-based filter via bcftools."""

    @dictify
    def _get_input_files_run(self, wildcards: Wildcards):
        yield from super()._get_input_files_run(wildcards).items()
        # The wrapper reads the BED file as input "include" or "exclude"
        yield self.config.regions.action, self.parent.get_upstream_paths("regions").bed

    def _get_params_run(self, wildcards: Wildcards) -> dict[str, Any]:
        return {"filter_name": self.parent.task_name, "mode": self.config.regions.mode}


class DkfzStepPart(_BamAwareStepPart):
    """DKFZ bias filter."""

    resource_usage = {
        "run": ResourceUsage(threads=1, runtime="12h", mem=f"{3 * 1024}MB"),
    }


class EbfilterStepPart(_BamAwareStepPart):
    """EBFilter."""

    actions = ("run", "write_panel")
    resource_usage = {
        "run": ResourceUsage(threads=1, runtime="24h", mem=f"{2 * 1024}MB"),
        "write_panel": ResourceUsage(threads=1, runtime="1h", mem=f"{2 * 1024}MB"),
    }

    def _get_input_files_write_panel(self, wildcards):
        return {}

    @dictify
    def _get_input_files_run(self, wildcards: Wildcards):
        yield from super()._get_input_files_run(wildcards).items()
        cfg: EbfilterConfig = self.config.ebfilter
        panel_file = cfg.path_panel_of_normals_sample_list
        if not panel_file:
            panel_file = self._get_output_files_write_panel()["txt"]
        yield "txt", panel_file

    def _get_output_files_write_panel(self):
        return {"txt": "work/eb_filter.panel_of_normals/out/eb_filter.panel_of_normals.txt"}

    def _get_params_run(self, wildcards: Wildcards) -> dict[str, Any]:
        params = super()._get_params_run(wildcards)
        params.pop("path_panel_of_normals_sample_list", None)
        return params | {
            "has_annotation": getattr(self.config, "has_annotation", True),
        }

    def write_panel_of_normals_file(self, wildcards):
        out_path = self._get_output_files_write_panel()["txt"]
        with open(out_path, "wt") as fh:
            for bam in self._get_panel_of_normal_bams(wildcards):
                print(bam, file=fh)

    @listify
    def _get_panel_of_normal_bams(self, wildcards):
        df = self.parent.build_library_dataframe()
        normal_df = df[(df["role"] == "normal") & (df["extraction_type"] == "dna")]
        libraries = sorted(normal_df["library_name"].tolist())

        cfg: EbfilterConfig = self.config.ebfilter
        random.seed(cfg.shuffle_seed)
        random.shuffle(libraries)
        for lib_name in libraries[: cfg.panel_of_normals_size]:
            aln: ExpectedAlignments = self.parent.get_upstream_paths(
                "alignments", library_name=lib_name
            )
            if isinstance(aln, dict):
                yield aln["bam"]
            else:
                yield aln.bam


# ---------------------------------------------------------------------------
# Map tool name → step part class
# ---------------------------------------------------------------------------

_TOOL_STEP_PART: dict[str, type[VariantFiltrationStepPart]] = {
    "bcftools": BcftoolsStepPart,
    "vembrane": VembraneStepPart,
    "regions": RegionsStepPart,
    "dkfz": DkfzStepPart,
    "ebfilter": EbfilterStepPart,
}


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

    @classmethod
    def task_produces(cls, config, upstream):
        """The variants this task reads, tagged ``filtered``."""
        return (upstream["variants"].with_tags("filtered"),)

    @classmethod
    def filter_ids(cls, config):
        """The IDs of the DKFZ bias filter, and the expression names of ``vembrane`` in tag mode."""
        if config.tool == "dkfz":
            return frozenset({"bPcr", "bSeq"})
        if config.tool == "vembrane" and config.vembrane.mode == "tag":
            return frozenset(config.vembrane.expressions)
        return frozenset()

    config_model_class = VariantFiltrationConfigModel

    default_relationships = {
        "matched_normal_lib": RelationshipDefinition(
            via="donor_name",
            target="role == 'normal' and extraction_type == 'dna'",
        )
    }

    @classmethod
    def get_output_paths(cls, config, signature=None, **kwargs) -> dict[str, str]:
        """Return local filtered-variant output paths for downstream consumers."""
        lib = kwargs.get("library_name", "{library_name}")
        return {
            "vcf": f"work/{lib}/out/{lib}.vcf.gz",
            "vcf_tbi": f"work/{lib}/out/{lib}.vcf.gz.tbi",
        }

    def __init__(self, workflow, project, task_name):
        super().__init__(workflow, project, task_name)
        tool_cls = _TOOL_STEP_PART[self.config.tool]
        self.register_sub_step_classes((tool_cls,))
