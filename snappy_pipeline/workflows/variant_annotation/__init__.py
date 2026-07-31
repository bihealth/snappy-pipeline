# -*- coding: utf-8 -*-
"""Implementation of the unified ``variant_annotation`` step.

One task = one annotation tool. Chain via ``depends_on.variant`` to annotate
somatic or germline VCF inputs with either VEP or Mehari.
"""

import os

from biomedsheets.shortcuts import GenericSampleSheet
from snakemake.io import expand

from snappy_pipeline.utils import dictify, listify
from snappy_pipeline.workflows.abstract import (
    BaseStep,
    BaseStepPart,
    LinkOutStepPart,
    ResourceUsage,
)
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType

from .model import VariantAnnotation as VariantAnnotationConfigModel, ExpectedVariantVcf

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

EXT_VALUES = (".vcf.gz", ".vcf.gz.tbi", ".vcf.gz.md5", ".vcf.gz.tbi.md5")
EXT_NAMES = ("vcf", "vcf_tbi", "vcf_md5", "vcf_tbi_md5")

_WORK_PREFIX = os.path.join("work", "{library_name}")
_OUT_PREFIX = os.path.join(_WORK_PREFIX, "out", "{library_name}")
_LOG_PREFIX = os.path.join(_WORK_PREFIX, "log", "{library_name}")


class VariantAnnotationStepPart(BaseStepPart):
    """Base class for concrete variant annotation tool step parts."""

    actions = ("run",)

    def get_input_files(self, action):
        self._validate_action(action)
        return self._get_input_files_run

    @dictify
    def _get_input_files_run(self, wildcards):
        lib = wildcards.library_name
        variant = ExpectedVariantVcf.model_validate(
            self.parent.get_upstream_paths("variant", library_name=lib)
        )
        yield "vcf", variant.vcf
        yield "vcf_tbi", variant.vcf_tbi
        yield "reference", self.w_config.static_data_config.reference.path

    @dictify
    def get_output_files(self, action):
        self._validate_action(action)
        for name, ext in zip(EXT_NAMES, EXT_VALUES):
            yield name, _OUT_PREFIX + ext

    @dictify
    def _get_log_file(self, action):
        self._validate_action(action)
        for key, ext in (
            ("log", ".log"),
            ("conda_info", ".conda_info.txt"),
            ("conda_list", ".conda_list.txt"),
        ):
            yield key, _LOG_PREFIX + ext

    def get_args(self, action):
        self._validate_action(action)
        return {"config": getattr(self.config, self.name).model_dump(by_alias=True)}


class VepStepPart(VariantAnnotationStepPart):
    name = "vep"

    def get_resource_usage(self, action: str, **kwargs) -> ResourceUsage:
        self._validate_action(action)
        num_threads = self.config.vep.num_threads
        return ResourceUsage(
            threads=num_threads,
            runtime="1d",
            mem=f"{2 * num_threads}GB",
        )


class MehariStepPart(VariantAnnotationStepPart):
    name = "mehari"

    @dictify
    def _get_input_files_run(self, wildcards):
        yield from super()._get_input_files_run(wildcards).items()
        if self.config.mehari and self.config.mehari.transcripts:
            yield "transcripts", self.config.mehari.transcripts
        if self.config.mehari and self.config.mehari.frequencies:
            yield "frequencies", self.config.mehari.frequencies
        if self.config.mehari and self.config.mehari.clinvar:
            yield "clinvar", self.config.mehari.clinvar

    def get_resource_usage(self, action: str, **kwargs) -> ResourceUsage:
        self._validate_action(action)
        return ResourceUsage(
            threads=self.config.mehari.threads,
            runtime="20m",
            mem="8GB",
        )


_TOOL_STEP_PART = {
    "vep": VepStepPart,
    "mehari": MehariStepPart,
}


class VariantAnnotationWorkflow(BaseStep):
    """Annotate germline or somatic variant calls with a single selected tool."""

    name = "variant_annotation"
    consumes = {DataSignature(DataType.VARIANTS): True}
    produces = [
        DataSignature(DataType.VARIANTS, frozenset({"germline", "snv", "indel", "annotated"})),
        DataSignature(DataType.VARIANTS, frozenset({"somatic", "snv", "indel", "annotated"})),
    ]

    config_model_class = VariantAnnotationConfigModel
    sheet_shortcut_class = GenericSampleSheet

    @classmethod
    def get_output_paths(cls, signature=None, **kwargs) -> dict[str, str]:
        """Return local annotated VCF output paths for somatic or germline variants."""
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
        for entity_name in self.output_entities:
            yield from expand(
                os.path.join("output", "{library_name}", "out", "{library_name}{ext}"),
                library_name=[entity_name],
                ext=EXT_VALUES,
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

    def check_config(self):
        self.ensure_w_config(
            ("static_data_config", "reference", "path"),
            "Path to reference FASTA not configured but required for variant annotation",
        )
