# -*- coding: utf-8 -*-
"""Implementation of the ``variant_filtration`` step

This step takes annotated variants as the input from ``variant_annotation`` and performs various
filtration and postprocessing operations:

1. filter to high-confidence variants
    1. apply quality filter sets
    2. filter for consistency between different callers
2. filter to compatible mode of inheritance
3. filter by population/cohort frequency, remove polymorphisms
4. filter by region
5. filter by scores (e.g., conservation)
6. filter for het. comp. inheritance or keep all

# ::

#     1
#     stringent
#     loose

#     2
#     $qual.denovo
#     $qual.dom
#     $qual.rec_hom

#     3
#     $qual.denovo.denov_freq
#     $qual.dom.dom_freq
#     $qual.dom.rec_freq
#     $qual.rec_hom.rec_freq

#     4
#     $qual.denovo.denov_freq.$region
#     $qual.dom.dom_freq.$region
#     $qual.dom.rec_freq.$region
#     $qual.rec_hom.rec_freq.$region

#     5
#     $qual.denovo.denov_freq.$region.$scores
#     $qual.dom.dom_freq.$region.$scores
#     $qual.dom.rec_freq.$region.$scores
#     $qual.rec_hom.rec_freq.$region.$scores

#     6
#     $qual.denovo.denov_freq.$region.keep_all
#     $qual.dom.dom_freq.$region.keep_all
#     $qual.dom.rec_freq.$region.$scores.same_gene
#     $qual.dom.rec_freq.$region.$scores.same_tad
#     $qual.dom.rec_freq.$region.$scores.itv_500bp
#     $qual.rec_hom.rec_freq.$region.keep_all

================
Filtration Steps
================

The combinations of the filters is given in the configuration setting ``filter_combinations``
as dot-separated values, e.g., ``AA.BB.CC``.

==========
Step Input
==========

TODO

===========
Step Output
===========

TODO

====================
Global Configuration
====================

TODO

=====================
Default Configuration
=====================

The default configuration is as follows.

.. include:: DEFAULT_CONFIG_variant_filtration.rst

=======
Reports
=======

Currently, no reports are generated.
"""

# TODO: the implementation is super ugly and needs some refinement...

import os
import os.path
import sys
from typing import Any

from biomedsheets.shortcuts import GermlineCaseSheet, is_not_background
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

from .model import (
    Ebfilter as EbfilterConfig,
    VariantFiltration as VariantFiltrationConfigModel,
)

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

EXT_VALUES = (".vcf.gz", ".vcf.gz.tbi", ".vcf.gz.md5", ".vcf.gz.tbi.md5")
EXT_NAMES = ("vcf", "vcf_tbi", "vcf_md5", "vcf_tbi_md5")

DEFAULT_CONFIG = VariantFiltrationConfigModel.default_config_yaml_string()

# Path template helpers
_WORK_PREFIX = os.path.join("work", "{library_name}")
_OUT_PREFIX = os.path.join(_WORK_PREFIX, "out", "{library_name}")
_LOG_PREFIX = os.path.join(_WORK_PREFIX, "log", "{library_name}")


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
        lib = wildcards.library_name
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
        # Build tumor→normal mapping from cancer sheets (if present).
        self._tumor_to_normal: dict[str, str] = {}
        for info, raw_sheet in zip(self.parent.data_set_infos, self.parent.sheets):
            if info.is_background or info.sheet_type != "matched_cancer":
                continue
            try:
                csheet = CancerCaseSheet(
                    raw_sheet,
                    options=CancerCaseSheetOptions(
                        allow_missing_normal=True, allow_missing_tumor=False
                    ),
                )
                for pair in csheet.all_sample_pairs:
                    t_lib = pair.tumor_sample.dna_ngs_library
                    n_lib = pair.normal_sample.dna_ngs_library if pair.normal_sample else None
                    if t_lib and n_lib:
                        self._tumor_to_normal[t_lib.name] = n_lib.name
            except Exception as exc:
                print(
                    f"WARNING: could not build tumor/normal mapping: {exc}",
                    file=sys.stderr,
                )

    @dictify
    def _get_input_files_run(self, wildcards: Wildcards):
        yield from super()._get_input_files_run(wildcards).items()

        yield "reference", self.w_config.static_data_config.reference.path

        lib = wildcards.library_name
        tumor_aln: ExpectedAlignments = self.parent.get_upstream_paths(
            "ngs_mapping", library_name=lib
        )
        if isinstance(tumor_aln, dict):
            yield "bam", tumor_aln["bam"]
        else:
            yield "bam", tumor_aln.bam

        normal_lib = self._tumor_to_normal.get(lib)
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
        libraries = []
        for info, raw_sheet in zip(self.parent.data_set_infos, self.parent.sheets):
            if info.sheet_type != "matched_cancer":
                continue
            try:
                csheet = CancerCaseSheet(
                    raw_sheet,
                    options=CancerCaseSheetOptions(
                        allow_missing_normal=True, allow_missing_tumor=False
                    ),
                )
                for donor in csheet.donors:
                    for bio_sample in donor.bio_samples.values():
                        if not bio_sample.extra_infos.get("isTumor", True):
                            if bio_sample.dna_ngs_library:
                                libraries.append(bio_sample.dna_ngs_library.name)
            except Exception:
                pass

        libraries.sort()
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

    @classmethod
    def default_config_yaml(cls):
        return DEFAULT_CONFIG

    @classmethod
    def get_output_paths(cls, signature=None, **kwargs) -> dict[str, str]:
        """Return local filtered-variant output paths for downstream consumers."""
        cls.require_signature(signature)
        lib = kwargs.get("library_name", "{library_name}")
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
        for lib_name in self._iter_library_names():
            yield from expand(
                os.path.join("output", "{library_name}", "out", "{library_name}{ext}"),
                library_name=[lib_name],
                ext=EXT_VALUES,
            )

    # ------------------------------------------------------------------
    # Library-name enumeration (sheet-type-aware)
    # ------------------------------------------------------------------

    @listify
    def _iter_library_names(self):
        """Yield relevant DNA library names from all non-background data sets."""
        for info, raw_sheet, shortcut_sheet in zip(
            self.data_set_infos, self.sheets, self.shortcut_sheets
        ):
            if info.is_background:
                continue
            if info.sheet_type == "matched_cancer":
                yield from self._cancer_library_names(raw_sheet)
            else:
                yield from self._generic_library_names(shortcut_sheet)

    @staticmethod
    def _cancer_library_names(raw_sheet):
        """Yield tumor DNA library names from a matched-cancer sheet."""
        try:
            csheet = CancerCaseSheet(
                raw_sheet,
                options=CancerCaseSheetOptions(
                    allow_missing_normal=True, allow_missing_tumor=False
                ),
            )
            for donor in csheet.donors:
                for bio_sample in donor.bio_samples.values():
                    if not bio_sample.extra_infos.get("isTumor", False):
                        continue
                    for ts in bio_sample.test_samples.values():
                        if ts.extra_infos.get("extractionType", "").lower() == "dna":
                            for lib in ts.ngs_libraries.values():
                                yield lib.name
        except Exception as exc:
            print(
                f"WARNING: could not enumerate cancer library names: {exc}",
                file=sys.stderr,
            )

    @staticmethod
    def _generic_library_names(shortcut_sheet):
        """Yield all DNA library names from a generic/germline sheet."""
        for lib in shortcut_sheet.all_ngs_libraries:
            ext = lib.test_sample.extra_infos.get("extractionType", "DNA")
            if ext.lower() == "dna":
                yield lib.name
