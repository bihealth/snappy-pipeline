# -*- coding: utf-8 -*-
"""Combine germline & somatic variants (useful for many applications)

==========
Step Input
==========

One ``vcf`` file for somatic variants, and another for germline variants.

Both these files can be annotated and/or filtered. But unfortunaterly, because steps
cannot appear multiple times in the pipeline, if annotation or filtration is used
for either somatic or germline file, then the combined vcf cannot be further annotated
or filtered (depending on which one has occurred).

===========
Step Output
===========

The combined ``vcf`` is found in ``output/{mapper}.combined.{tumor_library}/out/output/{mapper}.combined.{tumor_library}.vcf.gz``.
The variant callers, annotator(s) and filtration status are discarded in the filenames.
It means that for downstream applications (TMB, ...), the status of file must be un-filtered and
un-annotated, in order to generate the correct filenames.

====================
Global Configuration
====================

TODO

=====================
Default Configuration
=====================

The default configuration is as follows.

.. include:: DEFAULT_CONFIG_combine_variants.rst

=======
Reports
=======

Currently, no reports are generated.
"""

from typing import Any

from snakemake.io import Wildcards

from biomedsheets.shortcuts import CancerCaseSheet

from snappy_pipeline.base import MissingConfiguration
from snappy_pipeline.utils import dictify, listify
from snappy_pipeline.workflows.common.samplesheet import (
    filter_table_by_modality,
    sample_sheets,
    tumor_to_normal_mapping,
)
from snappy_pipeline.workflows.abstract import (
    BaseStep,
    BaseStepPart,
    LinkOutStepPart,
    ResourceUsage,
)

from snappy_pipeline.workflows.somatic_variant_calling import SomaticVariantCallingWorkflow
from snappy_pipeline.workflows.somatic_variant_annotation import SomaticVariantAnnotationWorkflow
from snappy_pipeline.workflows.somatic_variant_filtration import SomaticVariantFiltrationWorkflow

from snappy_pipeline.workflows.germline_variant_calling import GermlineVariantCallingWorkflow
from snappy_pipeline.workflows.germline_variant_annotation import GermlineVariantAnnotationWorkflow
from snappy_pipeline.workflows.germline_variant_filtration import GermlineVariantFiltrationWorkflow

from .model import CombineVariants as CombineVariantsConfigModel
from .model import InputVariantType
from .model import RenameCombine

__author__ = "Eric Blanc <eric.blanc@bih-charite.de>"

#: Default configuration for the any_variant_calling step
DEFAULT_CONFIG = CombineVariantsConfigModel.default_config_yaml_string()


class CombineVariantsStepPart(BaseStepPart):
    name = "combine"

    actions = ("run",)

    default_resource_usage = ResourceUsage(threads=1, memory="4G", time="03:59:59")

    def __init__(self, parent):
        super().__init__(parent)
        cfg: CombineVariantsConfigModel = self.config

        self.somatic_tpl = f"{cfg.tool_ngs_mapping}.{cfg.tool_somatic_variant_calling}"
        if annotator := cfg.tool_somatic_variant_annotation:
            self.somatic_tpl += f".{annotator}"
        if cfg.is_somatic_variant_filtered:
            self.somatic_tpl += ".filtered"
        self.somatic_tpl += ".{tumor_library}"

        self.germline_tpl = f"{cfg.tool_ngs_mapping}.{cfg.tool_germline_variant_calling}"
        if annotator := cfg.tool_germline_variant_annotation:
            self.germline_tpl += f".{annotator}"
        if cfg.is_germline_variant_filtered:
            self.germline_tpl += ".filtered"
        self.germline_tpl += ".{normal_library}"

        self.base_out = f"{cfg.tool_ngs_mapping}.combined.{{tumor_library}}"

    def get_input_files(self, action: str):
        self._validate_action(action)
        return getattr(self, f"_get_input_files_{action}")

    @dictify
    def _get_input_files_run(self, wildcards: Wildcards):
        yield "reference", self.w_config.static_data_config.reference.path

        vcf = "output/{tpl}/out/{tpl}.vcf.gz".format(tpl=self.somatic_tpl)
        somatic_variant = self.parent.sub_workflows["somatic_variant"]
        yield "somatic_vcf", somatic_variant(vcf.format(tumor_library=wildcards.tumor_library))

        vcf = "output/{tpl}/out/{tpl}.vcf.gz".format(tpl=self.germline_tpl)
        germline_variant = self.parent.sub_workflows["germline_variant"]
        yield (
            "germline_vcf",
            germline_variant(
                vcf.format(
                    normal_library=self.parent.tumor_to_normal_mapping[wildcards.tumor_library]
                )
            ),
        )

    def get_output_files(self, action: str) -> dict[str, Any]:
        match action:
            case "run":
                return {
                    k: "work/{tpl}/out/{tpl}.vcf.gz{ext}".format(tpl=self.base_out, ext=ext)
                    for k, ext in (
                        ("vcf", ""),
                        ("vcf_tbi", ".tbi"),
                        ("vcf_md5", ".md5"),
                        ("vcf_tbi_md5", ".tbi.md5"),
                    )
                }
            case _:
                raise MissingConfiguration(f"Unimplemented action {action}")

    def get_args(self, action: str):
        self._validate_action(action)
        return getattr(self, f"_get_args_{action}")

    def _get_args_run(self, wildcards: Wildcards) -> dict[str, Any]:
        if name_origin := self.config.rename_combined:
            if name_origin == RenameCombine.TUMOR:
                sample_name = wildcards.tumor_library
            elif name_origin == RenameCombine.GERMLINE:
                sample_name = self.parent.tumor_to_normal_mapping[wildcards.tumor_library]
            else:
                raise MissingConfiguration(f"Unimplemented sample name type {name_origin}")
            return {"tumor_library": wildcards.tumor_library, "sample_name": sample_name}
        else:
            return {"tumor_library": wildcards.tumor_library}

    @dictify
    def get_log_file(self, action: str):
        match action:
            case "run":
                for k, ext in (
                    ("log", "log"),
                    ("conda_list", "conda_list.txt"),
                    ("conda_info", "conda_info.txt"),
                ):
                    yield k, "work/{tpl}/log/{tpl}.{ext}".format(tpl=self.base_out, ext=ext)
                    yield (
                        k + "_md5",
                        "work/{tpl}/log/{tpl}.{ext}.md5".format(tpl=self.base_out, ext=ext),
                    )
            case _:
                raise MissingConfiguration(f"Unimplemented action {action}")


class CombineVariantsWorkflow(BaseStep):
    """Perform somatic variant filtration"""

    #: Workflow name
    name = "combine_variants"

    #: Default biomed sheet class
    sheet_shortcut_class = CancerCaseSheet

    @classmethod
    def default_config_yaml(cls):
        """Return default config YAML, to be overwritten by project-specific one."""
        return DEFAULT_CONFIG

    def __init__(self, workflow, config, config_lookup_paths, config_paths, workdir):
        previous_steps = []

        match config["step_config"][self.name]["somatic_variant_type"]:
            case InputVariantType.CALLING:
                previous_steps.append(SomaticVariantCallingWorkflow)
            case InputVariantType.ANNOTATION:
                previous_steps.append(SomaticVariantAnnotationWorkflow)
            case InputVariantType.FILTRATION:
                previous_steps.append(SomaticVariantFiltrationWorkflow)

        match config["step_config"][self.name]["germline_variant_type"]:
            case InputVariantType.CALLING:
                previous_steps.append(GermlineVariantCallingWorkflow)
            case InputVariantType.ANNOTATION:
                previous_steps.append(GermlineVariantAnnotationWorkflow)
            case InputVariantType.FILTRATION:
                previous_steps.append(GermlineVariantFiltrationWorkflow)

        super().__init__(
            workflow,
            config,
            config_lookup_paths,
            config_paths,
            workdir,
            config_model_class=CombineVariantsConfigModel,
            previous_steps=previous_steps,
        )
        # Register sub step classes so the sub steps are available
        self.register_sub_step_classes((CombineVariantsStepPart, LinkOutStepPart))

        # Register sub-workflows by variant origin
        self.register_sub_workflow(
            f"somatic_variant_{self.config['somatic_variant_type']}",
            self.config.path_somatic_variant,
            "somatic_variant",
        )
        self.register_sub_workflow(
            f"germline_variant_{self.config['germline_variant_type']}",
            self.config.path_germline_variant,
            "germline_variant",
        )

        self.table = filter_table_by_modality(sample_sheets(self.sheets), modality="dna")
        assert self.table.shape[1] > 0, "No valid samples"
        self.tumor_to_normal_mapping = tumor_to_normal_mapping(self.table)

    @listify
    def get_result_files(self):
        mapper = self.config.tool_ngs_mapping
        for t, n in self.tumor_to_normal_mapping.items():
            if n:
                for ext in ("", ".tbi", ".md5", ".tbi.md5"):
                    yield f"output/{mapper}.combined.{t}/out/{mapper}.combined.{t}.vcf.gz{ext}"
                for ext in ("log", "conda_list.txt", "conda_info.txt"):
                    for hash in ("", ".md5"):
                        yield f"output/{mapper}.combined.{t}/log/{mapper}.combined.{t}.{ext}{hash}"
