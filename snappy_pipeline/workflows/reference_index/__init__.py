# -*- coding: utf-8 -*-
"""Implementation of the ``reference_index`` step."""

from biomedsheets.shortcuts import GenericSampleSheet

from snappy_pipeline.utils import dictify, listify
from snappy_pipeline.workflows.abstract import BaseStep, BaseStepPart, ResourceUsage
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType

from .model import ReferenceIndex as ReferenceIndexConfigModel
from .model import Tool

DEFAULT_CONFIG = ReferenceIndexConfigModel.default_config_yaml_string()


class BuildReferenceCommonStepPart(BaseStepPart):
    name = "common"
    actions = ("run",)

    def get_resource_usage(self, action: str, **kwargs) -> ResourceUsage:
        self._validate_action(action)
        return ResourceUsage(threads=1, runtime="2h", mem="2GB")

    def get_input_files(self, action):
        self._validate_action(action)
        return {"reference": self.parent.get_reference_path()}

    @dictify
    def get_output_files(self, action):
        self._validate_action(action)
        work_prefix = "work/reference_index/out/reference"
        outputs = {
            "reference_fai": work_prefix + ".fa.fai",
            "reference_dict": work_prefix + ".dict",
            "reference_genome": work_prefix + ".fa.genome",
        }
        log_files = self.get_log_file(action)
        output_links = [p.replace("work/", "output/", 1) for p in outputs.values()]
        output_links.extend(p.replace("work/", "output/", 1) for p in log_files.values())
        yield from outputs.items()
        yield "output_links", output_links

    @dictify
    def _get_log_file(self, action):
        self._validate_action(action)
        prefix = "work/reference_index/log/reference_index.common"
        for key, ext in (
            ("log", ".log"),
            ("conda_info", ".conda_info.txt"),
            ("conda_list", ".conda_list.txt"),
        ):
            yield key, prefix + ext

    @dictify
    def get_args(self, action):
        self._validate_action(action)
        yield "reference", self.parent.get_reference_path()


class _IndexToolStepPart(BaseStepPart):
    output_suffixes: tuple[str, ...] = ()

    def get_input_files(self, action):
        self._validate_action(action)
        return {
            "reference": self.parent.get_reference_path(),
            "reference_fai": "work/reference_index/out/reference.fa.fai",
            "reference_dict": "work/reference_index/out/reference.dict",
            "reference_genome": "work/reference_index/out/reference.fa.genome",
        }

    @dictify
    def get_output_files(self, action):
        self._validate_action(action)
        work_prefix = "work/reference_index/out/reference"
        outputs = {}
        for suffix in self.output_suffixes:
            key = suffix.strip(".").replace(".", "_").replace("/", "_")
            outputs[key] = work_prefix + suffix

        log_files = self.get_log_file(action)
        output_links = [p.replace("work/", "output/", 1) for p in outputs.values()]
        output_links.extend(p.replace("work/", "output/", 1) for p in log_files.values())
        yield from outputs.items()
        yield "output_links", output_links

    @dictify
    def _get_log_file(self, action):
        self._validate_action(action)
        prefix = f"work/reference_index/log/reference_index.{self.name}"
        for key, ext in (
            ("log", ".log"),
            ("conda_info", ".conda_info.txt"),
            ("conda_list", ".conda_list.txt"),
        ):
            yield key, prefix + ext

    @dictify
    def get_args(self, action):
        self._validate_action(action)
        yield "reference", self.parent.get_reference_path()


class BwaIndexStepPart(_IndexToolStepPart):
    name = "bwa"
    actions = ("run",)
    output_suffixes = (".amb", ".ann", ".bwt", ".pac", ".sa")

    def get_resource_usage(self, action: str, **kwargs) -> ResourceUsage:
        self._validate_action(action)
        return ResourceUsage(threads=8, runtime="24h", mem="16GB")

    @dictify
    def get_args(self, action):
        yield from super().get_args(action).items()
        yield "algorithm", self.config.bwa.algorithm


class BwaMem2IndexStepPart(_IndexToolStepPart):
    name = "bwa_mem2"
    actions = ("run",)
    output_suffixes = (".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac")

    def get_resource_usage(self, action: str, **kwargs) -> ResourceUsage:
        self._validate_action(action)
        return ResourceUsage(threads=8, runtime="24h", mem="16GB")

    @dictify
    def get_args(self, action):
        yield from super().get_args(action).items()
        yield "extra_args", " ".join(self.config.bwa_mem2.extra_args)


class Minimap2IndexStepPart(_IndexToolStepPart):
    name = "minimap2"
    actions = ("run",)
    output_suffixes = (".mmi",)

    def get_resource_usage(self, action: str, **kwargs) -> ResourceUsage:
        self._validate_action(action)
        return ResourceUsage(threads=8, runtime="12h", mem="16GB")

    @dictify
    def get_args(self, action):
        yield from super().get_args(action).items()
        yield "extra_args", " ".join(self.config.minimap2.extra_args)


class StarIndexStepPart(_IndexToolStepPart):
    name = "star"
    actions = ("run",)
    output_suffixes = (".star/.done",)

    def get_resource_usage(self, action: str, **kwargs) -> ResourceUsage:
        self._validate_action(action)
        return ResourceUsage(threads=16, runtime="24h", mem="64GB")

    @dictify
    def get_args(self, action):
        yield from super().get_args(action).items()
        yield "extra_args", " ".join(self.config.star.extra_args)
        features = getattr(self.w_config.static_data_config, "features", None)
        yield "features", getattr(features, "path", "") if features else ""


class ReferenceIndexWorkflow(BaseStep):
    name = "reference_index"
    consumes = {}
    produces = [
        DataSignature(DataType.INDEX, frozenset({"bwa", "dna"})),
        DataSignature(DataType.INDEX, frozenset({"bwa_mem2", "dna"})),
        DataSignature(DataType.INDEX, frozenset({"minimap2", "dna"})),
        DataSignature(DataType.INDEX, frozenset({"star", "rna"})),
    ]

    sheet_shortcut_class = GenericSampleSheet
    config_model_class = ReferenceIndexConfigModel

    @classmethod
    def default_config_yaml(cls):
        return DEFAULT_CONFIG

    @classmethod
    def get_output_paths(cls, signature=None, **kwargs) -> dict[str, str]:
        cls.require_signature(signature)
        prefix = kwargs.get("prefix", "output/reference_index/out/reference")
        return {
            "bwa_index_prefix": prefix,
            "bwa_mem2_index_prefix": prefix,
            "minimap2_index": prefix + ".mmi",
            "star_index_dir": prefix + ".star",
            "reference_fai": prefix + ".fa.fai",
            "reference_dict": prefix + ".dict",
            "reference_genome": prefix + ".fa.genome",
        }

    def __init__(
        self,
        workflow,
        config,
        config_lookup_paths,
        config_paths,
        workdir,
        task_name: str = "",
        **kwargs,
    ):
        super().__init__(
            workflow,
            config,
            config_lookup_paths,
            config_paths,
            workdir,
            task_name=task_name,
            **kwargs,
        )
        tool_to_class = {
            Tool.bwa: BwaIndexStepPart,
            Tool.bwa_mem2: BwaMem2IndexStepPart,
            Tool.minimap2: Minimap2IndexStepPart,
            Tool.star: StarIndexStepPart,
        }
        selected_tool_class = tool_to_class[self.config.tool]
        self.register_sub_step_classes((BuildReferenceCommonStepPart, selected_tool_class))

    def get_reference_path(self) -> str:
        dep_task = getattr(self.config.depends_on, "reference_download", "")
        if dep_task:
            upstream = self.get_upstream_paths("reference_download")
            if getattr(upstream, "fasta", ""):
                return upstream.fasta
        return self.config.path_reference or self.w_config.static_data_config.reference.path

    @listify
    def get_result_files(self):
        yield from self.get_output_files("common", "run")["output_links"]
        yield from self.get_output_files(str(self.config.tool), "run")["output_links"]
