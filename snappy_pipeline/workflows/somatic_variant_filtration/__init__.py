# -*- coding: utf-8 -*-
"""Implementation of the ``somatic_variant_filtration`` step"""

import os
import random
import sys
from collections import OrderedDict
from typing import Any

from biomedsheets.shortcuts import CancerCaseSheet, CancerCaseSheetOptions, is_not_background
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
from snappy_pipeline.workflows.somatic_variant_calling.model import ExpectedSomaticVariants

from .model import Ebfilter as EbfilterConfig
from .model import SomaticVariantFiltration as SomaticVariantFiltrationConfigModel

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

EXT_VALUES = (".vcf.gz", ".vcf.gz.tbi", ".vcf.gz.md5", ".vcf.gz.tbi.md5")
EXT_NAMES = ("vcf", "vcf_tbi", "vcf_md5", "vcf_tbi_md5")

DEFAULT_CONFIG = SomaticVariantFiltrationConfigModel.default_config_yaml_string()


class SomaticVariantFiltrationStepPart(BaseStepPart):
    def __init__(self, parent):
        super().__init__(parent)
        self.config = parent.config
        self.name_pattern = "{tumor_library}"
        self.tumor_ngs_library_to_sample_pair = OrderedDict()
        for sheet in self.parent.shortcut_sheets:
            self.tumor_ngs_library_to_sample_pair.update(
                sheet.all_sample_pairs_by_tumor_dna_ngs_library
            )
        self.donors = OrderedDict()
        for sheet in self.parent.shortcut_sheets:
            for donor in sheet.donors:
                self.donors[donor.name] = donor
        self.tumor_to_normal_library = OrderedDict()
        for tumor_library, normal_sample in self.tumor_ngs_library_to_sample_pair.items():
            for test_sample in normal_sample.normal_sample.bio_sample.test_samples.values():
                for ngs_library in test_sample.ngs_libraries.values():
                    if tumor_library not in self.tumor_to_normal_library:
                        self.tumor_to_normal_library[tumor_library] = (
                            test_sample.name + "-" + ngs_library.secondary_id
                        )

    def get_normal_lib_name(self, wildcards):
        pair = self.tumor_ngs_library_to_sample_pair.get(wildcards.tumor_library, None)
        return pair.normal_sample.dna_ngs_library.name if pair else None


class OneFilterStepPart(SomaticVariantFiltrationStepPart):
    name = "one_filter"
    actions = ("run",)
    resource_usage = {"run": ResourceUsage(threads=1, runtime="2h", mem=f"{8 * 1024}MB")}

    def get_input_files(self, action):
        self._validate_action(action)
        return getattr(self, f"_get_input_files_{action}")

    @dictify
    def _get_input_files_run(self, wildcards, is_unpack=True):
        filter_nb = int(wildcards["filter_nb"])
        name_pattern = self.name_pattern.format(**wildcards)
        if filter_nb > 1:
            prev = list(self.config.filter_list[filter_nb - 2].keys())[0]
            n = filter_nb - 1
            yield (
                "vcf",
                os.path.join("work", name_pattern, "out", name_pattern + f".{prev}_{n}.vcf.gz"),
            )
        else:
            somatic_variant: ExpectedSomaticVariants = self.parent.get_upstream_paths(
                "somatic_variant", library_name=name_pattern
            )
            yield "vcf", somatic_variant.vcf

    @dictify
    def get_output_files(self, action):
        self._validate_action(action)
        prefix = os.path.join(
            "work",
            self.name_pattern,
            "out",
            self.name_pattern + "." + self.filter_name + "_{filter_nb}",
        )
        key_ext = {
            "vcf": ".vcf.gz",
            "vcf_tbi": ".vcf.gz.tbi",
            "vcf_md5": ".vcf.gz.md5",
            "vcf_tbi_md5": ".vcf.gz.tbi.md5",
        }
        for key, ext in key_ext.items():
            yield key, prefix + ext

    @dictify
    def get_log_file(self, action):
        self._validate_action(action)
        key_ext = (
            ("log", ".log"),
            ("conda_info", ".conda_info.txt"),
            ("conda_list", ".conda_list.txt"),
        )
        for key, ext in key_ext:
            yield (
                key,
                os.path.join(
                    "work",
                    self.name_pattern,
                    "log",
                    self.name_pattern + "." + self.filter_name + "_{filter_nb}" + ext,
                ),
            )
            yield (
                key + "_md5",
                os.path.join(
                    "work",
                    self.name_pattern,
                    "log",
                    self.name_pattern + "." + self.filter_name + "_{filter_nb}" + ext + ".md5",
                ),
            )

    def get_args(self, action):
        self._validate_action(action)
        return self._get_args

    def _get_args(self, wildcards: Wildcards) -> dict[str, Any]:
        filter_nb = int(wildcards["filter_nb"])
        params = self.config.filter_list[filter_nb - 1][self.filter_name].model_dump(by_alias=True)
        params["filter_name"] = "{}_{}".format(self.filter_name, wildcards["filter_nb"])
        return params


class OneFilterWithBamStepPart(OneFilterStepPart):
    @dictify
    def _get_input_files_run(self, wildcards, **_kwargs):
        parent = super(OneFilterWithBamStepPart, self)._get_input_files_run
        yield from parent(wildcards, **_kwargs).items()

        yield "reference", self.w_config.static_data_config.reference.path

        name_pattern = "{tumor_library}".format(**wildcards)
        tumor_alignments: ExpectedAlignments = self.parent.get_upstream_paths(
            "ngs_mapping", library_name=name_pattern
        )
        yield "bam", tumor_alignments.bam
        if normal_library := self.tumor_to_normal_library.get(wildcards["tumor_library"], None):
            normal_alignments: ExpectedAlignments = self.parent.get_upstream_paths(
                "ngs_mapping", library_name=normal_library
            )
            yield "normal", normal_alignments.bam


class OneFilterDkfzStepPart(OneFilterWithBamStepPart):
    name = "one_dkfz"
    filter_name = "dkfz"
    resource_usage = {"run": ResourceUsage(threads=1, runtime="12h", mem=f"{3 * 1024}MB")}


class OneFilterEbfilterStepPart(OneFilterWithBamStepPart):
    name = "one_ebfilter"
    filter_name = "ebfilter"
    actions = ("run", "write_panel")

    resource_usage = {
        "run": ResourceUsage(threads=1, runtime="24h", mem=f"{2 * 1024}MB"),
        "write_panel": ResourceUsage(threads=1, runtime="1h", mem=f"{2 * 1024}MB"),
    }

    @dictify
    def _get_input_files_run(self, wildcards, **_kwargs):
        parent = super(OneFilterEbfilterStepPart, self)._get_input_files_run
        yield from parent(wildcards, **_kwargs).items()
        cfg: EbfilterConfig = self._get_args(wildcards)
        sample_files = cfg["path_panel_of_normals_sample_list"]
        if not sample_files:
            sample_files = self._get_output_files_write_panel()["txt"].format(**wildcards)
        yield "txt", sample_files

    def _get_output_files_write_panel(self):
        return {"txt": "work/eb_filter.panel_of_normals/out/eb_filter.panel_of_normals.txt"}

    def get_output_files(self, action):
        output_files = super(OneFilterEbfilterStepPart, self).get_output_files(action)
        if action == "write_panel":
            output_files = self._get_output_files_write_panel()
        return output_files

    def _get_args(self, wildcards: Wildcards) -> dict[str, Any]:
        return super(OneFilterEbfilterStepPart, self)._get_args(wildcards) | {
            "has_annotation": self.config.has_annotation,
        }

    def write_panel_of_normals_file(self, wildcards):
        output_path = self.get_output_files("write_panel")["txt"].format(**wildcards)
        with open(output_path, "wt") as outf:
            for bam_path in self._get_panel_of_normal_bams(wildcards):
                print(bam_path, file=outf)

    @listify
    def _get_panel_of_normal_bams(self, wildcards):
        libraries = []
        for sheet in self.parent.shortcut_sheets:
            for donor in sheet.donors:
                for bio_sample in donor.bio_samples.values():
                    if not bio_sample.extra_infos["isTumor"]:
                        libraries.append(bio_sample.dna_ngs_library.name)
        libraries.sort()

        for filter_cfg in self.config.filter_list:
            if list(filter_cfg.keys())[0] == "ebfilter":
                cfg: EbfilterConfig = list(filter_cfg.values())[0]
                break
        random.seed(cfg.shuffle_seed)
        lib_count = cfg["panel_of_normals_size"]
        random.shuffle(libraries)
        for library in libraries[:lib_count]:
            alignments: ExpectedAlignments = self.parent.get_upstream_paths(
                "ngs_mapping", library_name=library
            )
            yield alignments.bam


class OneFilterBcftoolsStepPart(OneFilterStepPart):
    name = "one_bcftools"
    filter_name = "bcftools"


class OneFilterRegionsStepPart(OneFilterStepPart):
    name = "one_regions"
    filter_name = "regions"


class OneFilterVembraneStepPart(OneFilterStepPart):
    name = "one_vembrane"
    filter_name = "vembrane"


class OneFilterProtectedStepPart(OneFilterStepPart):
    name = "one_protected"
    filter_name = "protected"


class LastFilterStepPart(SomaticVariantFiltrationStepPart):
    name = "last_filter"
    actions = ("run",)

    def get_input_files(self, action):
        self._validate_action(action)

        filter_names = [list(filter_name.keys())[0] for filter_name in self.config.filter_list]
        filter_nb = len(self.config.filter_list)
        filter_name = filter_names[filter_nb - 1]
        vcf = os.path.join(
            "work",
            self.name_pattern,
            "out",
            self.name_pattern + f".{filter_name}_{filter_nb}.vcf.gz",
        )
        prefix = os.path.join("work", self.name_pattern, "log", self.name_pattern)
        logs = [
            prefix + "." + filter_name + "_" + str(filter_nb + 1) + "." + e + m
            for filter_nb, filter_name in enumerate(filter_names)
            for e in ("log", "conda_list.txt", "conda_info.txt")
            for m in ("", ".md5")
        ]
        return {"vcf": vcf, "logs": logs}

    @dictify
    def get_output_files(self, action):
        self._validate_action(action)
        name_pattern = "{tumor_library}"
        vcf = os.path.join("work", name_pattern, "out", name_pattern)
        merged_log = os.path.join("work", name_pattern, "log", name_pattern + ".merged.tar.gz")
        return {
            "vcf": vcf + ".vcf.gz",
            "vcf_tbi": vcf + ".vcf.gz.tbi",
            "vcf_md5": vcf + ".vcf.gz.md5",
            "vcf_tbi_md5": vcf + ".vcf.gz.tbi.md5",
            "full": vcf + ".full.vcf.gz",
            "full_tbi": vcf + ".full.vcf.gz.tbi",
            "full_md5": vcf + ".full.vcf.gz.md5",
            "full_tbi_md5": vcf + ".full.vcf.gz.tbi.md5",
            "log": merged_log,
            "log_md5": merged_log + ".md5",
        }

    @dictify
    def get_log_file(self, action):
        self._validate_action(action)
        name_pattern = "{tumor_library}"
        tpl = os.path.join("work", name_pattern, "log", name_pattern)
        return {
            "log": tpl + ".log",
            "log_md5": tpl + ".log.md5",
            "conda_list": tpl + ".conda_list.txt",
            "conda_list_md5": tpl + ".conda_list.txt.md5",
            "conda_info": tpl + ".conda_info.txt",
            "conda_info_md5": tpl + ".conda_info.txt.md5",
        }


class SomaticVariantFiltrationWorkflow(BaseStep):
    name = "somatic_variant_filtration"
    consumes = {
        DataSignature(
            DataType.VARIANTS, frozenset({"somatic", ("snv", "indel"), "annotated"})
        ): True
    }
    produces = [
        DataSignature(DataType.VARIANTS, frozenset({"somatic", "snv", "indel", "filtered"}))
    ]

    config_model_class = SomaticVariantFiltrationConfigModel
    sheet_shortcut_class = CancerCaseSheet

    sheet_shortcut_kwargs = {
        "options": CancerCaseSheetOptions(allow_missing_normal=True, allow_missing_tumor=False)
    }

    @classmethod
    def get_output_paths(cls, signature=None, **kwargs) -> dict[str, str]:
        """Return local filtered-variant output paths for downstream consumers."""
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
            (
                OneFilterDkfzStepPart,
                OneFilterEbfilterStepPart,
                OneFilterBcftoolsStepPart,
                OneFilterRegionsStepPart,
                OneFilterVembraneStepPart,
                OneFilterProtectedStepPart,
                LastFilterStepPart,
                LinkOutStepPart,
            )
        )
        # Inputs are resolved via get_upstream_paths() in step parts.

    @listify
    def get_result_files(self):
        if not self.config.filter_list:
            return  # nothing to filter → no output files
        log_ext = [e + m for e in ("log", "conda_list.txt", "conda_info.txt") for m in ("", ".md5")]
        name_pattern = "{tumor_library}"

        yield from self._yield_result_files_matched(
            os.path.join("output", name_pattern, "out", name_pattern + "{ext}"),
            ext=[f + e for f in ("", ".full") for e in EXT_VALUES],
        )
        yield from self._yield_result_files_matched(
            os.path.join("output", name_pattern, "log", name_pattern + ".{ext}"),
            ext=log_ext,
        )
        yield from self._yield_result_files_matched(
            os.path.join("output", name_pattern, "log", name_pattern + ".merged.tar.gz{ext}"),
            ext=("", ".md5"),
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
