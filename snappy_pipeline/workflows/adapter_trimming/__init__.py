# -*- coding: utf-8 -*-
"""Adapter trimming of FASTQ files

Trims adapters and low-quality bases from the reads of each library.
The ``tool`` field selects ``bbduk`` (BBDuk of BBTools) or ``fastp``.
Each library gives one result, which keeps the sub-directory and file names of its input FASTQ files.
Tasks that read reads name this task in ``depends_on.reads``.

Outputs
-------

Below ``tasks/<task>/output/``, the files of each library are in ``{library_name}/``:

- ``out/.done`` marks the trimmed FASTQ files, which are in the same directory.
- ``report/.done`` marks the tool reports (``bbduk``: statistics and histograms per input file, ``fastp``: HTML and JSON reports), which are in the same directory.

Reads that the tool filters out are kept below ``rejected/`` in ``work/`` and are not published.
Downstream tasks read the trimmed files from ``work/{library_name}/out/`` through the key ``fastq_dir``.
"""

import os

from biomedsheets.shortcuts import GenericSampleSheet

from snappy_pipeline.utils import dictify
from snappy_pipeline.workflows.abstract import (
    BaseStep,
    BaseStepPart,
    ResourceUsage,
)
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType
from snappy_pipeline.workflows.common.reads import reads_input_files, reads_params

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

    #: Keys of the tool config that name files, which are inputs and not params
    file_keys: tuple[str, ...] = ()

    def _get_input_files_run(self, wildcards):
        inputs = reads_input_files(self.parent, wildcards.library_name)
        tool_config = self.config.get(self.name)
        for key in self.file_keys:
            if value := getattr(tool_config, key):
                inputs[key] = value
        return inputs

    published_keys = ("out_done", "report_done")

    @dictify
    def _get_output_files_run(self):
        tool = self.config.tool
        if self.name != tool:
            return []
        return (
            ("out_done", "work/{library_name}/out/.done"),
            ("report_done", "work/{library_name}/report/.done"),
            ("rejected_done", "work/{library_name}/rejected/.done"),
        )

    @dictify
    def _get_log_file(self, action):
        tool = self.config.tool
        if self.name != tool:
            return []
        yield "log", "work/{library_name}/log/{library_name}.log"

    def _get_params_run(self, wildcards):
        # The trimmed files keep the sub-directory and name of each input file, in input order.
        output = {"left": [], "right": []}
        for group in self.parent.read_groups(wildcards.library_name):
            for mate, files in output.items():
                if mate in group.paths:
                    files.append(
                        {
                            "relative_path": os.path.dirname(group.relpaths[mate]) or ".",
                            "filename": os.path.basename(group.relpaths[mate]),
                        }
                    )
        return {
            "library_name": wildcards.library_name,
            "output_left": output["left"],
            "output_right": output["right"],
            "config": {
                key: value
                for key, value in dict(self.config.get(self.name)).items()
                if key not in self.file_keys
            },
        } | reads_params(self.parent, wildcards.library_name)


class BbdukStepPart(AdapterTrimmingStepPart):
    name = "bbduk"
    file_keys = ("adapter_sequences", "barcodes")

    @property
    def default_resource_usage(self) -> ResourceUsage:
        return ResourceUsage(threads=self.config.bbduk.num_threads, runtime="12h", mem="24000MB")


class FastpStepPart(AdapterTrimmingStepPart):
    name = "fastp"
    file_keys = ("filter_by_index1", "filter_by_index2")

    @property
    def default_resource_usage(self) -> ResourceUsage:
        return ResourceUsage(threads=self.config.fastp.num_threads, runtime="12h", mem="24000MB")


class AdapterTrimmingWorkflow(BaseStep):
    name = "adapter_trimming"
    produces = [DataSignature(DataType.RAW, frozenset({"trimmed"}))]

    sheet_shortcut_class = GenericSampleSheet
    config_model_class = AdapterTrimmingConfigModel

    @classmethod
    def get_output_paths(cls, config, signature=None, **kwargs) -> dict[str, str]:
        """Return local output paths for trimmed/raw FASTQ consumption."""
        return {"fastq_dir": "work"}

    def __init__(self, workflow, project, task_name):
        super().__init__(workflow, project, task_name)
        match self.config.tool:
            case "bbduk":
                selected = BbdukStepPart
            case "fastp":
                selected = FastpStepPart
            case _:
                raise NotImplementedError(f"Unknown tool: {self.config.tool}")
        self.register_sub_step_classes((selected,))
