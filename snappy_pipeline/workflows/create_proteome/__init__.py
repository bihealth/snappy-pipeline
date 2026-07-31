import os
from typing import Any

from biomedsheets.shortcuts import GermlineCaseSheet
from snakemake.io import expand
from snakemake.iocontainers import Wildcards

from snappy_pipeline.base import MissingConfiguration
from snappy_pipeline.utils import dictify, listify
from snappy_pipeline.workflows.abstract import (
    BaseStep,
    BaseStepPart,
    LinkOutStepPart,
    ResourceUsage,
)
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType
from snappy_pipeline.workflows.variant_annotation import VariantAnnotationWorkflow
from snappy_pipeline.workflows.variant_calling import VariantCallingWorkflow
from snappy_pipeline.workflows.variant_filtration import VariantFiltrationWorkflow

from .model import CreateProteome as CreateProteomeConfigModel

__author__ = "Eric Blanc <eric.blanc@bih-charite.de>"


_OUT_PREFIX = "work/{library_name}/out/{library_name}"
_LOG_PREFIX = "work/{library_name}/log/{library_name}"


class CreateProteomeStepPart(BaseStepPart):
    name = "create_proteome"
    actions = ("run",)
    default_resource_usage = ResourceUsage(threads=1, mem="4G", runtime="4h")

    def get_input_files(self, action: str):
        self._validate_action(action)
        return getattr(self, f"_get_input_files_{action}")

    @dictify
    def _get_input_files_run(self, wildcards: Wildcards):
        yield "reference", self.w_config.static_data_config.reference.path
        yield "features", self.w_config.static_data_config.features.path

        if self.config.path_proteome:
            yield "proteome", self.config.path_proteome

        variant = self.parent.get_upstream_paths("variant", library_name=wildcards.library_name)
        yield "vcf", getattr(variant, "vcf", None) or variant["vcf"]

    def get_output_files(self, action: str) -> dict[str, Any]:
        match action:
            case "run":
                return {"vcf": _OUT_PREFIX + ".fa.gz"}
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
                    ("log", ".log"),
                    ("conda_list", ".conda_list.txt"),
                    ("conda_info", ".conda_info.txt"),
                ):
                    yield k, _LOG_PREFIX + ext
                    yield k + "_md5", _LOG_PREFIX + ext + ".md5"
            case _:
                raise MissingConfiguration(f"Unimplemented action {action}")


class CreateProteomeWorkflow(BaseStep):
    name = "create_proteome"
    sheet_shortcut_class = GermlineCaseSheet

    consumes = {DataSignature(DataType.VARIANTS): True}
    produces = [DataSignature(DataType.TABULAR, frozenset({"proteome"}))]

    config_model_class = CreateProteomeConfigModel

    @classmethod
    def get_output_paths(cls, signature=None, **kwargs) -> dict[str, str]:
        cls.require_signature(signature)
        lib = kwargs.get("library_name", "{library_name}")
        return {
            "proteome": f"output/{lib}/out/{lib}.fa.gz",
        }

    def __init__(self, workflow, config, config_lookup_paths, config_paths, workdir, **kwargs):
        previous_steps = [
            VariantCallingWorkflow,
            VariantAnnotationWorkflow,
            VariantFiltrationWorkflow,
        ]

        super().__init__(
            workflow,
            config,
            config_lookup_paths,
            config_paths,
            workdir,
            previous_steps=previous_steps,
            **kwargs,
        )
        self.register_sub_step_classes((CreateProteomeStepPart, LinkOutStepPart))

    @listify
    def get_result_files(self):
        for entity_name in self.output_entities:
            yield from expand(
                os.path.join("output", "{library_name}", "out", "{library_name}.fa.gz{hash}"),
                library_name=[entity_name],
                hash=("", ".md5"),
            )
            yield from expand(
                os.path.join("output", "{library_name}", "log", "{library_name}{ext}"),
                library_name=[entity_name],
                ext=(
                    ".log",
                    ".log.md5",
                    ".conda_info.txt",
                    ".conda_info.txt.md5",
                    ".conda_list.txt",
                    ".conda_list.txt.md5",
                ),
            )
