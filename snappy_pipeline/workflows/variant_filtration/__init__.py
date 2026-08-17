# -*- coding: utf-8 -*-
"""Implementation of the unified ``variant_filtration`` step.

One task = one filter tool.  Chain multiple tasks via ``depends_on.variant``
to build a sequential filter pipeline.

Supported tools
---------------
bcftools    – expression-based filter tagging
vembrane    – expression-based tagging or filtering via mode switch
regions     – region/BED-based filter tagging
dkfz        – DKFZ bias filter (requires aligned BAMs, cancer sheet)
ebfilter    – EBFilter (requires aligned BAMs, cancer sheet)
"""

import os
import random
from typing import Any

from biomedsheets.shortcuts import GenericSampleSheet
from snakemake.io import expand
from snakemake.iocontainers import Wildcards

from snappy_pipeline.utils import dictify, listify
from snappy_pipeline.workflows.abstract import (
    BaseStep,
    BaseStepPart,
    LinkOutStepPart,
    ResourceUsage,
)
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType
from snappy_pipeline.workflows.ngs_mapping.model import ExpectedAlignments
from snappy_pipeline.models import RelationshipDefinition

from .model import (
    Ebfilter as EbfilterConfig,
    VariantFiltration as VariantFiltrationConfigModel,
)

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

EXT_VALUES = (".vcf.gz", ".vcf.gz.tbi", ".vcf.gz.md5", ".vcf.gz.tbi.md5")
EXT_NAMES = ("vcf", "vcf_tbi", "vcf_md5", "vcf_tbi_md5")


# Path template helpers
_WORK_PREFIX = os.path.join("work", "{tumor_library}")
_OUT_PREFIX = os.path.join(_WORK_PREFIX, "out", "{tumor_library}")
_LOG_PREFIX = os.path.join(_WORK_PREFIX, "log", "{tumor_library}")


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

    def get_input_files(self, action):
        self._validate_action(action)
        return self._get_input_files_run

    @dictify
    def _get_input_files_run(self, wildcards: Wildcards):
        lib = wildcards.tumor_library
        variant = self.parent.get_upstream_paths("variant", library_name=lib)
        # Accept both typed schema and plain dict
        if isinstance(variant, dict):
            yield "vcf", variant["vcf"]
            yield "vcf_tbi", variant["vcf_tbi"]
        else:
            yield "vcf", variant.vcf
            yield "vcf_tbi", variant.vcf_tbi

    @dictify
    def get_output_files(self, action):
        self._validate_action(action)
        for name, ext in zip(EXT_NAMES, EXT_VALUES):
            yield name, _OUT_PREFIX + ext

    @dictify
    def get_log_file(self, action):
        self._validate_action(action)
        for key, ext in (
            ("log", ".log"),
            ("log_md5", ".log.md5"),
            ("conda_info", ".conda_info.txt"),
            ("conda_info_md5", ".conda_info.txt.md5"),
            ("conda_list", ".conda_list.txt"),
            ("conda_list_md5", ".conda_list.txt.md5"),
        ):
            yield key, _LOG_PREFIX + ext

    def get_args(self, action):
        self._validate_action(action)
        return self._get_args

    def _get_args(self, wildcards: Wildcards) -> dict[str, Any]:
        cfg = getattr(self.config, self.config.tool)
        params: dict[str, Any] = {
            # Keep wrapper naming stable and task-scoped for chaining.
            "filter_name": getattr(self, "filter_name", self.parent.task_name or self.config.tool),
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

        yield "reference", self.w_config.static_data_config.reference.path

        lib = wildcards.tumor_library
        tumor_aln: ExpectedAlignments = self.parent.get_upstream_paths(
            "ngs_mapping", library_name=lib
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
                "ngs_mapping", library_name=normal_lib
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

    filter_name = "bcftools"


class VembraneStepPart(VariantFiltrationStepPart):
    """Unified vembrane tag/filter step."""

    filter_name = "vembrane"

    def _get_args(self, wildcards: Wildcards) -> dict[str, Any]:
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
                    "aux": cfg.aux,
                    "context": cfg.context,
                    "context_files": cfg.context_files,
                    "ontology": cfg.ontology,
                }
            )
        return params


class RegionsStepPart(VariantFiltrationStepPart):
    """Region/BED-based filter via bcftools."""

    filter_name = "regions"


class DkfzStepPart(_BamAwareStepPart):
    """DKFZ bias filter."""

    filter_name = "dkfz"

    resource_usage = {
        "run": ResourceUsage(threads=1, runtime="12h", mem=f"{3 * 1024}MB"),
    }


class EbfilterStepPart(_BamAwareStepPart):
    """EBFilter."""

    filter_name = "ebfilter"

    actions = ("run", "write_panel")
    resource_usage = {
        "run": ResourceUsage(threads=1, runtime="24h", mem=f"{2 * 1024}MB"),
        "write_panel": ResourceUsage(threads=1, runtime="1h", mem=f"{2 * 1024}MB"),
    }

    def get_input_files(self, action):
        self._validate_action(action)
        if action == "write_panel":
            return {}
        return self._get_input_files_run

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

    def get_output_files(self, action):
        self._validate_action(action)
        if action == "write_panel":
            return self._get_output_files_write_panel()
        return super().get_output_files(action)

    def _get_args(self, wildcards: Wildcards) -> dict[str, Any]:
        return super()._get_args(wildcards) | {
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
                "ngs_mapping", library_name=lib_name
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
    """Unified variant filtration – applies a single configurable filter tool to any VCF input.

    Supports somatic-specific tools (dkfz, ebfilter) as well as generic tools
    (bcftools, vembrane, regions). Use multiple tasks with
    ``depends_on.variant`` to compose a sequential filter pipeline.
    """

    name = "variant_filtration"
    consumes = {DataSignature(DataType.VARIANTS): True}
    produces = [DataSignature(DataType.VARIANTS, frozenset({"filtered"}))]

    config_model_class = VariantFiltrationConfigModel
    sheet_shortcut_class = GenericSampleSheet

    default_relationships = {
        "matched_normal_lib": RelationshipDefinition(
            via="donor_name",
            target="role == 'normal' and extraction_type == 'dna'",
        )
    }

    @classmethod
    def get_output_paths(cls, signature=None, **kwargs) -> dict[str, str]:
        """Return local filtered-variant output paths for downstream consumers."""
        cls.require_signature(signature)
        lib = kwargs.get("tumor_library", "{tumor_library}")
        return {
            "vcf": f"output/{lib}/out/{lib}.vcf.gz",
            "vcf_tbi": f"output/{lib}/out/{lib}.vcf.gz.tbi",
        }

    def __init__(
        self,
        workflow,
        config,
        config_lookup_paths,
        config_paths,
        workdir,
        task_name: str | None = None,
        **kwargs,
    ):
        super().__init__(
            workflow,
            config,
            config_lookup_paths,
            config_paths,
            workdir,
            previous_steps=(),
            task_name=task_name,
            **kwargs,
        )
        tool_cls = _TOOL_STEP_PART[self.config.tool]
        self.register_sub_step_classes((tool_cls, LinkOutStepPart))

    @listify
    def get_result_files(self):
        for entity_name in self.output_entities:
            yield from expand(
                os.path.join("output", "{tumor_library}", "out", "{tumor_library}{ext}"),
                tumor_library=[entity_name],
                ext=EXT_VALUES,
            )
