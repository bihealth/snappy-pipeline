# -*- coding: utf-8 -*-
"""Implementation of the ``adapter_trimming`` step"""

import os

from biomedsheets.shortcuts import GenericSampleSheet
from snakemake.io import expand

from snappy_pipeline.utils import dictify, listify
from snappy_pipeline.workflows.abstract import (
    BaseStep,
    BaseStepPart,
    ResourceUsage,
)
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType

from .model import AdapterTrimming as AdapterTrimmingConfigModel

#: Adatper trimming tools
TRIMMERS = ("bbduk", "fastp")

#: Default configuration for the hla_typing schema


class AdapterTrimmingStepPart(BaseStepPart):
    """Adapter trimming common features"""

    #: Step name
    name = ""

    #: Class available actions
    actions = ("run",)

    def __init__(self, parent):
        super().__init__(parent)
        self.base_path_out = "work/{library_name}"

    @dictify
    def _get_input_files_run(self, wildcards):
        yield "reads", self.parent.reads_input(wildcards.library_name)

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

    def _get_params_run(self, wildcards):
        # The trimmed files keep the sub-directory and name of each input file.
        reads = {"reads_left": {}, "reads_right": {}}
        for group in self.parent.read_groups(wildcards.library_name):
            for mate in ("left", "right"):
                if mate in group.paths:
                    reads[f"reads_{mate}"][group.paths[mate]] = {
                        "relative_path": os.path.dirname(group.relpaths[mate]) or ".",
                        "filename": os.path.basename(group.relpaths[mate]),
                    }
        return {
            "library_name": wildcards.library_name,
            "input": reads,
            "config": dict(self.config.get(self.name)),
        }


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
    produces = [DataSignature(DataType.RAW, frozenset({"trimmed"}))]

    sheet_shortcut_class = GenericSampleSheet
    config_model_class = AdapterTrimmingConfigModel

    @classmethod
    def get_output_paths(cls, config, signature=None, **kwargs) -> dict[str, str]:
        """Return local output paths for trimmed/raw FASTQ consumption."""
        return {"fastq_dir": "output"}

    def __init__(self, workflow, project, task_name):
        super().__init__(workflow, project, task_name)
        match self.config.tool:
            case "bbduk":
                selected = BbdukStepPart
            case "fastp":
                selected = FastpStepPart
            case _:
                raise NotImplementedError(f"Unknown tool: {self.config.tool}")
        self.register_sub_step_classes((LinkOutFastqStepPart, selected))

    @listify
    def get_result_files(self):
        tpls = (
            "output/{library_name}/out/.done",
            "output/{library_name}/report/.done",
            "output/{library_name}/log/.done",
        )
        for library_name in self.output_entities:
            for tpl in tpls:
                yield tpl.format(library_name=library_name)
