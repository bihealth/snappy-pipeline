# -*- coding: utf-8 -*-
"""Implementation of the ``somatic_variant_calling`` step"""

import os
import sys
from collections import OrderedDict
from itertools import chain

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

from .model import (
    SomaticVariantCalling as SomaticVariantCallingConfigModel,
)
from .model import TumorNormalMode as TumorNormalMode

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

EXT_VALUES = (".vcf.gz", ".vcf.gz.tbi", ".vcf.gz.md5", ".vcf.gz.tbi.md5")
EXT_NAMES = ("vcf", "vcf_tbi", "vcf_md5", "vcf_tbi_md5")

EXT_MATCHED = {
    "mutect": {
        "vcf": ".vcf.gz",
        "vcf_md5": ".vcf.gz.md5",
        "vcf_tbi": ".vcf.gz.tbi",
        "vcf_tbi_md5": ".vcf.gz.tbi.md5",
        "full_vcf": ".full.vcf.gz",
        "full_vcf_md5": ".full.vcf.gz.md5",
        "full_vcf_tbi": ".full.vcf.gz.tbi",
        "full_vcf_tbi_md5": ".full.vcf.gz.tbi.md5",
        "txt": ".txt",
        "txt_md5": ".txt.md5",
        "wig": ".wig",
        "wig_md5": ".wig.md5",
    },
    "scalpel": {
        "vcf": ".vcf.gz",
        "vcf_md5": ".vcf.gz.md5",
        "vcf_tbi": ".vcf.gz.tbi",
        "vcf_tbi_md5": ".vcf.gz.tbi.md5",
        "full_vcf": ".full.vcf.gz",
        "full_vcf_md5": ".full.vcf.gz.md5",
        "full_vcf_tbi": ".full.vcf.gz.tbi",
        "full_vcf_tbi_md5": ".full.vcf.gz.tbi.md5",
        "tar": ".tar.gz",
        "tar_md5": ".tar.gz.md5",
    },
    "mutect2": {
        "vcf": ".vcf.gz",
        "vcf_md5": ".vcf.gz.md5",
        "vcf_tbi": ".vcf.gz.tbi",
        "vcf_tbi_md5": ".vcf.gz.tbi.md5",
        "full_vcf": ".full.vcf.gz",
        "full_vcf_md5": ".full.vcf.gz.md5",
        "full_vcf_tbi": ".full.vcf.gz.tbi",
        "full_vcf_tbi_md5": ".full.vcf.gz.tbi.md5",
    },
}

SOMATIC_VARIANT_CALLERS = {"mutect2"}
DEFAULT_CONFIG = SomaticVariantCallingConfigModel.default_config_yaml_string()


class SomaticVariantCallingStepPart(BaseStepPart):
    def __init__(self, parent):
        super().__init__(parent)
        self.base_path_out = "work/{{tumor_library}}/out/{{tumor_library}}{ext}"
        self.tumor_ngs_library_to_sample_pair = OrderedDict()
        for sheet in self.parent.shortcut_sheets:
            self.tumor_ngs_library_to_sample_pair.update(
                sheet.all_sample_pairs_by_tumor_dna_ngs_library
            )

    def get_input_files(self, action: str):
        self._validate_action(action)
        return getattr(self, f"_get_input_files_{action}")

    @dictify
    def _get_input_files_run(self, wildcards: Wildcards):
        ngs_mapping = self.parent.modules["ngs_mapping"]
        tumor_base_path = ("output/{tumor_library}/out/{tumor_library}").format(**wildcards)

        input_files = {
            "tumor_bam": ngs_mapping(tumor_base_path + ".bam"),
            "tumor_bai": ngs_mapping(tumor_base_path + ".bam.bai"),
        }

        normal_library = self.get_normal_lib_name(wildcards)
        if normal_library:
            normal_base_path = "output/{normal_library}/out/{normal_library}".format(
                normal_library=normal_library, **wildcards
            )
            input_files.update(
                {
                    "normal_bam": ngs_mapping(normal_base_path + ".bam"),
                    "normal_bai": ngs_mapping(normal_base_path + ".bam.bai"),
                }
            )

        return input_files

    def get_normal_lib_name(self, wildcards):
        pair = self.tumor_ngs_library_to_sample_pair.get(wildcards.tumor_library, None)
        return pair.normal_sample.dna_ngs_library.name if pair else None

    def get_tumor_lib_name(self, wildcards):
        pair = self.tumor_ngs_library_to_sample_pair.get(wildcards.tumor_library, None)
        return pair.tumor_sample.dna_ngs_library.name if pair else wildcards.tumor_library

    def get_output_files(self, action):
        self._validate_action(action)
        return dict(zip(EXT_NAMES, expand(self.base_path_out, ext=EXT_VALUES)))

    @dictify
    def _get_log_file(self, action):
        self._validate_action(action)

        prefix = "work/{{tumor_library}}/log/{{tumor_library}}"
        key_ext = (
            ("log", ".log"),
            ("conda_info", ".conda_info.txt"),
            ("conda_list", ".conda_list.txt"),
        )
        for key, ext in key_ext:
            yield key, prefix + ext
            yield key + "_md5", prefix + ext + ".md5"


class Mutect2StepPart(SomaticVariantCallingStepPart):
    name = "mutect2"

    actions = [
        "scatter",
        "run",
        "gather",
        "filter",
    ]

    resource_usage_dict = {
        "scatter": ResourceUsage(threads=1, runtime="2m", mem="1000MB"),
        "run": ResourceUsage(threads=1, runtime="5d", mem="8000MB"),
        "gather": ResourceUsage(threads=1, runtime="4h", mem="32768MB"),
        "filter": ResourceUsage(threads=2, runtime="4h", mem="15872MB"),
        "contamination": ResourceUsage(threads=2, runtime="4h", mem="7680MB"),
        "pileup_normal": ResourceUsage(threads=2, runtime="4h", mem="8000MB"),
        "pileup_tumor": ResourceUsage(threads=2, runtime="4h", mem="8000MB"),
    }

    def __init__(self, parent):
        super().__init__(parent)
        run_resource_usage = self.resource_usage_dict["run"]
        self.resource_usage_dict["run"] = ResourceUsage(
            threads=self.config.mutect2.num_cores or run_resource_usage.threads,
            runtime=run_resource_usage.runtime,
            mem=run_resource_usage.mem,
        )

    def check_config(self):
        tool = self.config.tool
        if self.name != tool:
            return
        self.parent.ensure_w_config(
            ("static_data_config", "reference", "path"),
            "Path to reference FASTA not configured but required for %s" % (self.name,),
        )

    def get_input_files(self, action):
        self._validate_action(action)
        return getattr(self, "_get_input_files_{}".format(action))

    def get_args(self, action):
        self._validate_action(action)
        return getattr(self, f"_get_args_{action}")

    def _get_args_scatter(self, wildcards):
        ignore_chroms = list(
            set(
                self.w_config.get("ignore_chroms", [])
                + self.config.get("ignore_chroms", [])
                + self.config.get(self.name).get("ignore_chroms", [])
            )
        )
        return {
            "ignore_chroms": sorted(list(ignore_chroms)),
            "padding": self.config.mutect2.padding,
            "java_options": self.config.mutect2.contamination.java_options,
            "extra_arguments": self.config.mutect2.contamination.extra_arguments,
        }

    def _get_args_pileup_normal(self, wildcards):
        return self.config.mutect2.contamination.pileup.model_dump(by_alias=True) | {
            "normal_lib_name": self.get_normal_lib_name(wildcards)
        }

    def _get_args_pileup_tumor(self, wildcards):
        return self.config.mutect2.contamination.pileup.model_dump(by_alias=True) | {
            "tumor_lib_name": self.get_tumor_lib_name(wildcards)
        }

    def _get_args_contamination(self, wildcards):
        return {
            "java_options": self.config.mutect2.contamination.java_options,
            "extra_arguments": self.config.mutect2.contamination.extra_arguments,
        }

    def _get_args_run(self, wildcards):
        return {
            "normal_lib_name": self.get_normal_lib_name(wildcards),
            "java_options": self.config.mutect2.java_options,
            "extra_arguments": self.config.mutect2.extra_arguments,
        }

    def _get_args_gather(self, wildcards):
        return {}

    def _get_args_filter(self, wildcards):
        return self.config.mutect2.filtration.model_dump(by_alias=True)

    def _get_input_files_scatter(self, wildcards):
        return {"fai": self.w_config.static_data_config.reference.path + ".fai"}

    def _get_input_files_run(self, wildcards):
        ngs_mapping = self.parent.modules["ngs_mapping"]
        tumor_base_path = ("output/{tumor_library}/out/{tumor_library}").format(**wildcards)
        scatteritem_base_path = (
            "work/{tumor_library}/out/{tumor_library}/mutect2par/scatter/{scatteritem}".format(
                **wildcards
            )
        )

        input_files = {
            "tumor_bam": ngs_mapping(tumor_base_path + ".bam"),
            "tumor_bai": ngs_mapping(tumor_base_path + ".bam.bai"),
            "region": scatteritem_base_path + ".region.bed",
        }

        tumor_normal_mode = self.config.mutect2.tumor_normal_mode
        if tumor_normal_mode != TumorNormalMode.TUMOR_ONLY:
            normal_library = self.get_normal_lib_name(wildcards)
            if normal_library:
                normal_base_path = "output/{normal_library}/out/{normal_library}".format(
                    normal_library=normal_library, **wildcards
                )
                input_files.update(
                    {
                        "normal_bam": ngs_mapping(normal_base_path + ".bam"),
                        "normal_bai": ngs_mapping(normal_base_path + ".bam.bai"),
                    }
                )
            else:
                if tumor_normal_mode == TumorNormalMode.PAIRED:
                    raise ValueError(
                        f"Normal sample for tumor {wildcards.tumor_library} required but not found."
                    )

        input_files["reference"] = self.w_config.static_data_config.reference.path

        if self.config.mutect2.germline_resource:
            input_files["germline_resource"] = self.config.mutect2.germline_resource
        if self.config.mutect2.panel_of_normals:
            input_files["panel_of_normals"] = self.config.mutect2.panel_of_normals

        return input_files

    def _get_input_files_gather(self, wildcards):
        gather = self.parent.workflow.globals.get("gather")
        gather = getattr(gather, self.name)
        scatteritem_base_path = (
            "work/{tumor_library}/out/{tumor_library}/mutect2par/run/{{scatteritem}}".format(
                **wildcards
            )
        )
        input_files = {
            "vcf": scatteritem_base_path + ".raw.vcf.gz",
            "stats": scatteritem_base_path + ".raw.vcf.stats",
            "f1r2": scatteritem_base_path + ".raw.f1r2.tar.gz",
        }
        return dict(map(lambda item: (item[0], gather(item[1])), input_files.items()))

    def _get_input_files_filter(self, wildcards):
        base_path = "work/{tumor_library}/out/{tumor_library}".format(**wildcards)
        input_files = {
            "raw": base_path + ".raw.vcf.gz",
            "stats": base_path + ".raw.vcf.stats",
            "orientation": base_path + ".raw.read_orientation_model.tar.gz",
            "reference": self.w_config.static_data_config.reference.path,
        }
        if self.get_normal_lib_name(wildcards):
            if self.config.mutect2.contamination.enabled:
                input_files["table"] = base_path + ".contamination.tbl"
                input_files["segments"] = base_path + ".segments.tbl"
        return input_files

    def _get_input_files_pileup_normal(self, wildcards):
        ngs_mapping = self.parent.modules["ngs_mapping"]
        base_path = "output/{normal_library}/out/{normal_library}".format(
            normal_library=self.get_normal_lib_name(wildcards), **wildcards
        )
        return {
            "bam": ngs_mapping(base_path + ".bam"),
            "bai": ngs_mapping(base_path + ".bam"),
            "reference": self.w_config.static_data_config.reference.path,
            "common_variants": self.config.mutect2.contamination.common_variants,
        }

    def _get_input_files_pileup_tumor(self, wildcards):
        ngs_mapping = self.parent.modules["ngs_mapping"]
        base_path = "output/{tumor_library}/out/{tumor_library}".format(**wildcards)
        return {
            "bam": ngs_mapping(base_path + ".bam"),
            "bai": ngs_mapping(base_path + ".bam"),
            "reference": self.w_config.static_data_config.reference.path,
            "common_variants": self.config.mutect2.contamination.common_variants,
        }

    def _get_input_files_contamination(self, wildcards: Wildcards):
        base_path = "work/{tumor_library}/out/{tumor_library}".format(**wildcards)
        return {
            "normal": base_path + ".normal.pileup",
            "tumor": base_path + ".tumor.pileup",
            "reference": self.w_config.static_data_config.reference.path,
        }

    def get_output_files(self, action):
        exts = {}
        output_files = {}

        self._validate_action(action)
        tool = self.config.tool
        if self.name != tool:
            return {}
        base_path_out = self.base_path_out

        if action == "scatter":
            scatter = self.parent.workflow.globals.get("scatter")
            scatter = getattr(scatter, self.name)
            template = "work/{{tumor_library}}/out/{{tumor_library}}/mutect2par/scatter/{scatteritem}.region.bed"
            return {"regions": scatter(template)}

        if action == "run":
            base_path_out = (
                "work/{{tumor_library}}/out/{{tumor_library}}/mutect2par/run/{{scatteritem}}{ext}"
            )
            exts = {
                "vcf": ".raw.vcf.gz",
                "vcf_md5": ".raw.vcf.gz.md5",
                "vcf_tbi": ".raw.vcf.gz.tbi",
                "vcf_tbi_md5": ".raw.vcf.gz.tbi.md5",
                "stats": ".raw.vcf.stats",
                "stats_md5": ".raw.vcf.stats.md5",
                "f1r2": ".raw.f1r2.tar.gz",
                "f1r2_md5": ".raw.f1r2.tar.gz.md5",
            }
        if action == "gather":
            exts = {
                "vcf": ".raw.vcf.gz",
                "vcf_md5": ".raw.vcf.gz.md5",
                "vcf_tbi": ".raw.vcf.gz.tbi",
                "vcf_tbi_md5": ".raw.vcf.gz.tbi.md5",
                "stats": ".raw.vcf.stats",
                "stats_md5": ".raw.vcf.stats.md5",
                "orientation": ".raw.read_orientation_model.tar.gz",
                "orientation_md5": ".raw.read_orientation_model.tar.gz.md5",
            }
        if action == "filter":
            exts = {
                "full_vcf": ".full.vcf.gz",
                "full_vcf_md5": ".full.vcf.gz.md5",
                "full_vcf_tbi": ".full.vcf.gz.tbi",
                "full_vcf_tbi_md5": ".full.vcf.gz.tbi.md5",
                "vcf": ".vcf.gz",
                "vcf_md5": ".vcf.gz.md5",
                "vcf_tbi": ".vcf.gz.tbi",
                "vcf_tbi_md5": ".vcf.gz.tbi.md5",
            }
        if action == "contamination":
            exts = {
                "table": ".contamination.tbl",
                "table_md5": ".contamination.tbl.md5",
                "segments": ".segments.tbl",
                "segments_md5": ".segments.tbl.md5",
            }
        if action == "pileup_normal":
            exts = {"pileup": ".normal.pileup", "pileup_md5": ".normal.pileup.md5"}
        if action == "pileup_tumor":
            exts = {"pileup": ".tumor.pileup", "pileup_md5": ".tumor.pileup.md5"}

        for k, v in exts.items():
            output_files[k] = base_path_out.format(ext=v)
        return output_files

    def get_log_file(self, action):
        postfix = ""
        log_files = {}
        key_ext = (
            ("log", ".log"),
            ("conda_info", ".conda_info.txt"),
            ("conda_list", ".conda_list.txt"),
        )

        self._validate_action(action)
        tool = self.config.tool
        if self.name != tool:
            return {}

        if action != "gather":
            if action == "run":
                postfix = ".{{scatteritem}}"
            else:
                postfix = "." + action

        prefix = ("work/{{tumor_library}}/log/{{tumor_library}}{postfix}").format(postfix=postfix)

        for key, ext in key_ext:
            log_files[key] = prefix + ext
            log_files[key + "_md5"] = prefix + ext + ".md5"
        return log_files

    def get_resource_usage(self, action: str, **kwargs) -> ResourceUsage:
        self._validate_action(action)
        return self.resource_usage_dict.get(action)


class SomaticVariantCallingWorkflow(BaseStep):
    name = "somatic_variant_calling"
    consumes = {DataSignature(DataType.ALIGNMENTS, frozenset({"dna"})): True}
    produces = [DataSignature(DataType.VARIANTS, frozenset({"somatic", "snv", "indel"}))]

    config_model_class = SomaticVariantCallingConfigModel
    sheet_shortcut_class = CancerCaseSheet
    sheet_shortcut_kwargs = {
        "options": CancerCaseSheetOptions(allow_missing_normal=True, allow_missing_tumor=True)
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
        tool = self.config.tool
        self.register_sub_step_classes(
            (
                Mutect2StepPart,
                LinkOutStepPart,
            )
        )
        self.register_module("ngs_mapping", "ngs_mapping")
        if tool == "mutect2":
            if self.config.mutect2.contamination.enabled:
                actions = self.sub_steps["mutect2"].actions
                self.sub_steps["mutect2"].actions = tuple(
                    chain(actions, ["contamination", "pileup_normal", "pileup_tumor"])
                )

    @listify
    def get_result_files(self):
        name_pattern = "{tumor_library.name}"
        caller = self.config.tool
        if caller in SOMATIC_VARIANT_CALLERS:
            yield from self._yield_result_files_matched(
                os.path.join("output", name_pattern, "out", name_pattern + "{ext}"),
                ext=EXT_MATCHED[caller].values() if caller in EXT_MATCHED else EXT_VALUES,
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
                            yield from expand(tpl, tumor_library=[ngs_library], **kwargs)
