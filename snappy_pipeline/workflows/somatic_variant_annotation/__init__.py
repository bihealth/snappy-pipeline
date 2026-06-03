# -*- coding: utf-8 -*-
"""Implementation of the ``somatic_variant_annotation`` step"""

import os
import sys
from collections import OrderedDict
from typing import cast

from biomedsheets.shortcuts import CancerCaseSheet, CancerCaseSheetOptions, is_not_background
from snakemake.io import expand

from snappy_pipeline.utils import dictify, listify
from snappy_pipeline.workflows.abstract import BaseStep, BaseStepPart, LinkOutStepPart
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType
from snappy_pipeline.workflows.ngs_mapping import ResourceUsage

from .model import (
    ExpectedVariantVcf,
    SomaticVariantAnnotation as SomaticVariantAnnotationConfigModel,
)

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

EXT_VALUES = (".vcf.gz", ".vcf.gz.tbi", ".vcf.gz.md5", ".vcf.gz.tbi.md5")
EXT_NAMES = ("vcf", "vcf_tbi", "vcf_md5", "vcf_tbi_md5")

ANNOTATION_TOOLS = ("vep", "mehari")
DEFAULT_CONFIG = SomaticVariantAnnotationConfigModel.default_config_yaml_string()


class AnnotateSomaticVcfStepPart(BaseStepPart):
    has_full = False

    def __init__(self, parent):
        super().__init__(parent)
        self.tumor_ngs_library_to_sample_pair = OrderedDict()
        for sheet in self.parent.shortcut_sheets:
            self.tumor_ngs_library_to_sample_pair.update(
                sheet.all_sample_pairs_by_tumor_dna_ngs_library
            )
        self.donors = OrderedDict()
        for sheet in self.parent.shortcut_sheets:
            for donor in sheet.donors:
                self.donors[donor.name] = donor

    @staticmethod
    def _name_template(
        config: SomaticVariantAnnotationConfigModel, annotator: str | None = None
    ) -> str:
        return "{tumor_library}"

    @dictify
    def get_input_files(self, action):
        self._validate_action(action)
        lib = self._name_template(self.config)
        variants = cast(
            ExpectedVariantVcf,
            self.parent.get_upstream_paths("variant", library_name=lib),
        )
        yield "vcf", variants.vcf
        yield "vcf_tbi", variants.vcf_tbi

    @dictify
    def get_output_files(self, action):
        self._validate_action(action)
        if self.name != self.config.tool:
            return []
        tpl = self._name_template(self.config)
        prefix = os.path.join("work", tpl, "out", tpl)
        key_ext = {"vcf": ".vcf.gz", "vcf_tbi": ".vcf.gz.tbi"}
        if self.has_full:
            key_ext["full"] = ".full.vcf.gz"
            key_ext["full_tbi"] = ".full.vcf.gz.tbi"
        for key, ext in key_ext.items():
            yield key, prefix + ext
            yield key + "_md5", prefix + ext + ".md5"

    @dictify
    def _get_log_file(self, action):
        self._validate_action(action)
        if self.name != self.config.tool:
            return []
        tpl = self._name_template(self.config, annotator=self.annotator)
        prefix = os.path.join("work", tpl, "log", tpl)

        key_ext = (
            ("log", ".log"),
            ("conda_info", ".conda_info.txt"),
            ("conda_list", ".conda_list.txt"),
        )
        for key, ext in key_ext:
            yield key, prefix + ext


class VepAnnotateSomaticVcfStepPart(AnnotateSomaticVcfStepPart):
    name = "vep"
    annotator = "vep"
    actions = ("run",)
    has_full = True

    PICK_ORDER = (
        "biotype",
        "mane_select",
        "mane_plus_clinical",
        "appris",
        "tsl",
        "ccds",
        "canonical",
        "rank",
        "length",
    )

    @dictify
    def get_input_files(self, action: str):
        input_files = super().get_input_files(action)
        for k, v in input_files.items():
            yield k, v
        yield "reference", self.w_config.static_data_config.reference.path

    def get_args(self, action):
        self._validate_action(action)

        def args_function(wildcards):
            return {"config": self.config.get(self.name).model_dump(by_alias=True)}

        return args_function

    def get_resource_usage(self, action: str, **kwargs) -> ResourceUsage:
        self._validate_action(action)
        return ResourceUsage(
            threads=self.config.vep.num_threads,
            runtime="24h",
            mem=f"{16 * 1024 * 1}MB",
        )


class MehariAnnotateSomaticVcfStepPart(AnnotateSomaticVcfStepPart):
    name = "mehari"
    annotator = "mehari"
    actions = ("run",)

    @dictify
    def get_input_files(self, action: str):
        input_files = super().get_input_files(action)
        for k, v in input_files.items():
            yield k, v

        yield "reference", self.w_config.static_data_config.reference.path

        if self.config.mehari and self.config.mehari.transcripts:
            yield "transcripts", self.config.mehari.transcripts
        if self.config.mehari and self.config.mehari.frequencies:
            yield "frequencies", self.config.mehari.frequencies
        if self.config.mehari and self.config.mehari.clinvar:
            yield "clinvar", self.config.mehari.clinvar

    def get_args(self, action):
        self._validate_action(action)

        def args_function(wildcards):
            return {"config": self.config.get(self.name).model_dump(by_alias=True)}

        return args_function

    def get_resource_usage(self, action: str, **kwargs) -> ResourceUsage:
        self._validate_action(action)
        return ResourceUsage(
            threads=self.config.mehari.threads,
            runtime="20m",
            mem="8GB",
        )


class SomaticVariantAnnotationWorkflow(BaseStep):
    name = "somatic_variant_annotation"
    consumes = {DataSignature(DataType.VARIANTS): True}
    produces = [
        DataSignature(DataType.VARIANTS, frozenset({"somatic", "snv", "indel", "annotated"}))
    ]

    config_model_class = SomaticVariantAnnotationConfigModel

    sheet_shortcut_class = CancerCaseSheet
    sheet_shortcut_kwargs = {
        "options": CancerCaseSheetOptions(allow_missing_normal=True, allow_missing_tumor=True)
    }

    @classmethod
    def get_output_paths(cls, signature=None, **kwargs) -> dict[str, str]:
        """Return local annotated VCF output paths for a somatic-variants signature."""
        cls.require_signature(signature)
        lib = kwargs.get("library_name", "{library_name}")
        return {
            "vcf": f"output/{lib}/out/{lib}.vcf.gz",
            "vcf_tbi": f"output/{lib}/out/{lib}.vcf.gz.tbi",
        }

    @classmethod
    def default_config_yaml(cls):
        return DEFAULT_CONFIG

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
        self.register_sub_step_classes(
            (VepAnnotateSomaticVcfStepPart, MehariAnnotateSomaticVcfStepPart, LinkOutStepPart)
        )

    @listify
    def get_result_files(self):
        if str(self.config.tool) not in ANNOTATION_TOOLS:
            return

        name_pattern = AnnotateSomaticVcfStepPart._name_template(self.config)
        yield from self._yield_result_files_matched(
            os.path.join("output", name_pattern, "out", name_pattern + "{ext}"),
            ext=EXT_VALUES,
        )
        yield from self._yield_result_files_matched(
            os.path.join("output", name_pattern, "log", name_pattern + "{ext}"),
            ext=(
                ".log",
                ".log.md5",
                ".conda_info.txt",
                ".conda_info.txt.md5",
                ".conda_list.txt",
                ".conda_list.txt.md5",
            ),
        )

        if self.sub_steps[self.config.tool].has_full:
            yield from self._yield_result_files_matched(
                os.path.join("output", name_pattern, "out", name_pattern + ".full{ext}"),
                ext=EXT_VALUES,
            )

    def _yield_result_files_matched(self, tpl, **kwargs):
        for sheet in filter(is_not_background, self.shortcut_sheets):
            for bio_entity in sheet.sheet.bio_entities.values():
                for bio_sample in bio_entity.bio_samples.values():
                    if not bio_sample.extra_infos.get("isTumor", False):
                        continue
                    for test_sample in bio_sample.test_samples.values():
                        extraction_type = test_sample.extra_infos.get("extractionType", "unknown")
                        if extraction_type.lower() != "dna":
                            if extraction_type == "unknown":
                                msg = "INFO: sample {} has missing extraction type, ignored"
                                print(msg.format(test_sample.name), file=sys.stderr)
                            continue
                        for ngs_library in test_sample.ngs_libraries.values():
                            yield from expand(tpl, tumor_library=[ngs_library.name], **kwargs)

    def _yield_result_files_joint(self, tpl, **kwargs):
        for sheet in filter(is_not_background, self.shortcut_sheets):
            for donor in sheet.donors:
                yield from expand(tpl, tumor_library=[donor], **kwargs)
