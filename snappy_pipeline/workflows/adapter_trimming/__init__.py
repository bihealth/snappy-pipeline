# -*- coding: utf-8 -*-
"""Implementation of the ``adapter_trimming`` step"""

import os
from collections import OrderedDict

from biomedsheets.shortcuts import GenericSampleSheet
from snakemake.io import expand

from snappy_pipeline.utils import dictify, listify
from snappy_pipeline.workflows.abstract import (
    BaseStep,
    BaseStepPart,
    LinkInPathGenerator,
    LinkInStepPart,
    ResourceUsage,
    get_ngs_library_folder_name,
)
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType

from .model import AdapterTrimming as AdapterTrimmingConfigModel

#: Adatper trimming tools
TRIMMERS = ("bbduk", "fastp")

#: Default configuration for the hla_typing schema
DEFAULT_CONFIG = AdapterTrimmingConfigModel.default_config_yaml_string()


class AdapterTrimmingStepPart(BaseStepPart):
    """Adapter trimming common features"""

    #: Step name
    name = ""

    #: Class available actions
    actions = ("run",)

    def __init__(self, parent):
        super().__init__(parent)
        self.base_path_in = "work/input_links/{library_name}"
        self.base_path_out = "work/{library_name}"
        #: Path generator for linking in
        self.path_gen = LinkInPathGenerator(
            self.parent.work_dir,
            self.parent.data_set_infos,
            self.parent.config_lookup_paths,
            preprocessed_path=self.parent.get_preprocessed_path(),
        )

    @dictify
    def get_input_files(self, action):
        self._validate_action(action)
        yield "done", "work/input_links/{library_name}/.done"

    @dictify
    def get_output_files(self, action):
        self._validate_action(action)
        tool = self.config.tool
        if self.name != tool:
            return []
        return (
            ("out_done", self.base_path_out + "/out/.done"),
            ("report_done", self.base_path_out + "/report/.done"),
            ("rejected_done", self.base_path_out + "/rejected/.done"),
        )

    @dictify
    def _get_log_file(self, action):
        self._validate_action(action)
        tool = self.config.tool
        if self.name != tool:
            return []
        _ = action
        prefix = "work/{library_name}/log/{library_name}"
        key_ext = (
            ("log", ".log"),
            ("conda_info", ".conda_info.txt"),
            ("conda_list", ".conda_list.txt"),
        )
        yield (
            "done",
            "work/{library_name}/log/.done",
        )
        for key, ext in key_ext:
            yield key, prefix + ext
            yield key + "_md5", prefix + ext + ".md5"

    def get_args(self, action):
        def args_function(wildcards):
            folder_name = get_ngs_library_folder_name(self.parent.sheets, wildcards.library_name)
            if self.parent.get_preprocessed_path():
                folder_name = wildcards.library_name
            reads_left = self._collect_reads(wildcards, folder_name, "")
            reads_right = self._collect_reads(wildcards, folder_name, "right-")
            return {
                "library_name": wildcards.library_name,
                "input": {
                    "reads_left": {key: reads_left[key] for key in sorted(reads_left.keys())},
                    "reads_right": {key: reads_right[key] for key in sorted(reads_right.keys())},
                },
                "config": dict(self.config.get(self.name)),
            }

        self._validate_action(action)
        return args_function

    def _collect_reads(self, wildcards, folder_name, prefix):
        task_prefix = self.parent.task_path_prefix()

        pattern_set_keys = ("right",) if prefix.startswith("right-") else ("left",)
        path_info = {}
        for _, path_infix, filename in self.path_gen.run(folder_name, pattern_set_keys):
            input_path = os.path.join(self.base_path_in, path_infix, filename).format(**wildcards)
            input_path = task_prefix + input_path

            assert input_path not in path_info.keys()
            paths = {
                "relative_path": path_infix,
                "filename": filename,
            }
            path_info[input_path] = paths
        return path_info


class BbdukStepPart(AdapterTrimmingStepPart):
    name = "bbduk"

    def get_resource_usage(self, action: str, **kwargs) -> ResourceUsage:
        self._validate_action(action)
        return ResourceUsage(
            threads=self.config.bbduk.num_threads,
            runtime="12h",
            mem="24000MB",
        )


class FastpStepPart(AdapterTrimmingStepPart):
    name = "fastp"

    def get_resource_usage(self, action: str, **kwargs) -> ResourceUsage:
        self._validate_action(action)
        return ResourceUsage(
            threads=self.config.fastp.num_threads,
            runtime="12h",
            mem="24000MB",
        )


class LinkOutFastqStepPart(BaseStepPart):
    name = "link_out_fastq"

    def __init__(self, parent):
        super().__init__(parent)
        self.base_path_in = "work/{wildcards.library_name}/{{sub_dir}}/.done"
        self.base_path_out = "output/{{library_name}}/{sub_dir}/.done"
        self.sub_dirs = ["log", "report", "out"]

    def get_input_files(self, action):
        def input_function(wildcards):
            return expand(self.base_path_in.format(wildcards=wildcards), sub_dir=self.sub_dirs)

        self._validate_action(action)
        return input_function

    def get_output_files(self, action):
        self._validate_action(action)
        return expand(self.base_path_out, sub_dir=self.sub_dirs)

    def run_locally(self, action, wildcards):
        self._validate_action(action)
        for sub_dir in self.sub_dirs:
            in_ = os.path.dirname(
                self.base_path_in.format(wildcards=wildcards).format(sub_dir=sub_dir)
            )
            out = os.path.dirname(self.base_path_out.format(sub_dir=sub_dir).format(**wildcards))
            os.makedirs(out, exist_ok=True)
            for root, d_names, f_names in os.walk(in_):
                rel_path = os.path.relpath(root, start=in_)
                for d_name in d_names:
                    os.makedirs(os.path.join(out, d_name), exist_ok=True)
                for f_name in f_names:
                    f = os.path.join(root, f_name)
                    if not os.path.islink(f):
                        target = os.path.relpath(f, start=os.path.join(out, rel_path))
                        os.symlink(target, os.path.join(out, rel_path, f_name))

    def _validate_action(self, action):
        assert action == "run"


class AdapterTrimmingWorkflow(BaseStep):
    name = "adapter_trimming"
    consumes = {DataSignature(DataType.RAW): True}
    produces = [DataSignature(DataType.RAW, frozenset({"trimmed"}))]

    sheet_shortcut_class = GenericSampleSheet
    config_model_class = AdapterTrimmingConfigModel

    @classmethod
    def get_output_paths(cls, signature=None, **kwargs) -> dict[str, str]:
        """Return local output paths for trimmed/raw FASTQ consumption."""
        cls.require_signature(signature)
        return {"fastq_dir": "output"}

    def __init__(self, *args, task_name: str, **kwargs):
        super().__init__(*args, task_name=task_name, **kwargs)
        self.register_sub_step_classes(
            (BbdukStepPart, FastpStepPart, LinkInStepPart, LinkOutFastqStepPart)
        )
        self.ngs_library_name_to_ngs_library = OrderedDict()
        for sheet in self.shortcut_sheets:
            for ngs_library in sheet.all_ngs_libraries:
                self.ngs_library_name_to_ngs_library[ngs_library.name] = ngs_library

    @classmethod
    def default_config_yaml(cls):
        return DEFAULT_CONFIG

    @listify
    def get_result_files(self):
        tpls = (
            "output/{ngs_library_name}/out/.done",
            "output/{ngs_library_name}/report/.done",
            "output/{ngs_library_name}/log/.done",
        )
        for sheet in self.shortcut_sheets:
            for ngs_library in sheet.all_ngs_libraries:
                for tpl in tpls:
                    yield tpl.format(ngs_library_name=ngs_library.name)
