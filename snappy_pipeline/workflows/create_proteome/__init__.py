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

The combined ``vcf`` is found in ``output/{mapper}.{caller}.{library}/out/output/{mapper}.combined.{library}.vcf.gz``.
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

.. include:: DEFAULT_CONFIG_create_proteome.rst

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
from snappy_pipeline.workflows.common.samplesheet import filter_table_by_modality, sample_sheets
from snappy_pipeline.workflows.abstract import (
    BaseStep,
    BaseStepPart,
    LinkOutStepPart,
    ResourceUsage,
)

from .model import CreateProteome as CreateProteomeConfigModel
from .model import InputVariantType

__author__ = "Eric Blanc <eric.blanc@bih-charite.de>"

#: Default configuration for the any_variant_calling step
DEFAULT_CONFIG = CreateProteomeConfigModel.default_config_yaml_string()


class CreateProteomeStepPart(BaseStepPart):
    name = "create_proteome"

    actions = ("run",)

    default_resource_usage = ResourceUsage(threads=1, memory="4G", time="03:59:59")

    def __init__(self, parent):
        super().__init__(parent)
        cfg: CreateProteomeConfigModel = self.config

        self.tpl = f"{cfg.tool_ngs_mapping}.{cfg.tool_variant_calling}"
        if annotator := cfg.tool_variant_annotation:
            self.tpl += f".{annotator}"
        if cfg.is_filtered:
            self.tpl += ".filtered"
        self.tpl += ".{library}"

    def get_input_files(self, action: str):
        self._validate_action(action)
        return getattr(self, f"_get_input_files_{action}")

    @dictify
    def _get_input_files_run(self, wildcards: Wildcards):
        yield "reference", self.w_config.static_data_config.reference.path
        yield "features", self.w_config.static_data_config.features.path

        if self.config.path_proteome:
            yield "proteome", self.config.path_proteome

        vcf = "output/{tpl}/out/{tpl}.vcf.gz".format(tpl=self.tpl)
        variant = self.parent.sub_workflows["variant"]
        yield "vcf", variant(vcf.format(library=wildcards.library))

    def get_output_files(self, action: str) -> dict[str, Any]:
        match action:
            case "run":
                return {"vcf": "work/{tpl}/out/{tpl}.vcf.gz".format(tpl=self.tpl)}
            case _:
                raise MissingConfiguration(f"Unimplemented action {action}")

    def get_args(self, action: str):
        self._validate_action(action)
        return getattr(self, f"_get_args_{action}")

    def _get_args_run(self, wildcards: Wildcards) -> dict[str, Any]:
        return {"add_reference": self.config.add_reference}

    @dictify
    def get_log_file(self, action: str):
        match action:
            case "run":
                for k, ext in (
                    ("log", "log"),
                    ("conda_list", "conda_list.txt"),
                    ("conda_info", "conda_info.txt"),
                ):
                    yield k, "work/{tpl}/log/{tpl}.{ext}".format(tpl=self.tpl, ext=ext)
                    yield (
                        k + "_md5",
                        "work/{tpl}/log/{tpl}.{ext}.md5".format(tpl=self.tpl, ext=ext),
                    )
            case _:
                raise MissingConfiguration(f"Unimplemented action {action}")


class CreateProteomeWorkflow(BaseStep):
    """Perform somatic variant filtration"""

    #: Workflow name
    name = "create_proteome"

    #: Default biomed sheet class
    sheet_shortcut_class = CancerCaseSheet

    @classmethod
    def default_config_yaml(cls):
        """Return default config YAML, to be overwritten by project-specific one."""
        return DEFAULT_CONFIG

    def __init__(self, workflow, config, config_lookup_paths, config_paths, workdir):
        previous_steps = []

        match config["step_config"][self.name]["variant_type"]:
            case InputVariantType.CALLING:
                from snappy_pipeline.workflows.any_variant_calling import AnyVariantCallingWorkflow

                previous_steps.append(AnyVariantCallingWorkflow)
            case InputVariantType.ANNOTATION:
                from snappy_pipeline.workflows.any_variant_annotation import (
                    AnyVariantAnnotationWorkflow,
                )

                previous_steps.append(AnyVariantAnnotationWorkflow)
            case InputVariantType.FILTRATION:
                from snappy_pipeline.workflows.any_variant_filtration import (
                    AnyVariantFiltrationWorkflow,
                )

                previous_steps.append(AnyVariantFiltrationWorkflow)

        super().__init__(
            workflow,
            config,
            config_lookup_paths,
            config_paths,
            workdir,
            config_model_class=CreateProteomeConfigModel,
            previous_steps=previous_steps,
        )
        # Register sub step classes so the sub steps are available
        self.register_sub_step_classes((CreateProteomeStepPart, LinkOutStepPart))

        # Register sub-workflows by variant origin
        self.register_sub_workflow(
            f"any_variant_{self.config['variant_type']}",
            self.config.path_variant,
            "variant",
        )

        self.table = filter_table_by_modality(sample_sheets(self.sheets), modality="dna")
        assert self.table.shape[1] > 0, "No valid samples"

    @listify
    def get_result_files(self):
        tpl = f"{self.config['tool_ngs_mapping']}.{self.config['tool_variant_calling']}"
        if self.config["tool_variant_annotation"]:
            tpl += f".{self.config['tool_variant_annotation']}"
        if self.config["is_filtered"]:
            tpl += ".filtered"
        for library in self.table["ngs_library"]:
            for hash in ("", ".md5"):
                yield f"output/{tpl}.{library}/out/{tpl}.{library}.fa.gz{hash}"
                for ext in ("log", "conda_list.txt", "conda_info.txt"):
                    yield f"output/{tpl}.{library}/log/{tpl}.{library}.{ext}{hash}"
