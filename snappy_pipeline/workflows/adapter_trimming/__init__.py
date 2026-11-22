# -*- coding: utf-8 -*-
"""Adapter trimming of FASTQ files

Trims adapters and low-quality bases from the reads of each library.
The ``tool`` field selects ``bbduk`` (BBDuk of BBTools) or ``fastp``.
``bbduk`` trims the sequences of the files in ``adapter_sequences`` and of the BBDuk built-in references in ``builtin_references``, such as ``adapters`` and ``phix``.
``bbduk`` needs at least one of the two fields.
Each library gives one result, whose trimmed files keep the sub-directory and file names of the input FASTQ files.
A task that reads the trimmed reads names this task in ``depends_on.reads``.

Outputs
-------

For each library ``{library_name}``, the task links the following result directories below ``tasks/<task>/output/`` as symlinks to ``work/``:

- ``{library_name}/out/`` holds the trimmed FASTQ files.
- ``{library_name}/report/`` holds the tool reports (``bbduk``: statistics and histograms per input file, ``fastp``: HTML and JSON reports).

Reads that the tool filters out stay in ``work/{library_name}/rejected/`` and are not result files.
Downstream tasks read the trimmed files from ``work/{library_name}/out/`` through the key ``fastq_dir``.
"""

import os

from snakemake.io import directory

from snappy_pipeline.utils import dictify
from snappy_pipeline.workflows.abstract import (
    BaseStep,
    BaseStepPart,
    ReportOutput,
    ResourceUsage,
)
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType
from snappy_pipeline.workflows.common.reads import reads_input_files, reads_params

from .model import AdapterTrimming as AdapterTrimmingConfigModel
from .model import Tool


class AdapterTrimmingStepPart(BaseStepPart):
    """Adapter trimming common features"""

    #: Step name
    name = ""

    #: Class available actions
    actions = ("run",)

    result_keys = ("out", "report")

    #: Keys of the tool config that name files, which are inputs and not params
    file_keys: tuple[str, ...] = ()

    def _get_log_file(self, action):
        return {"log": "work/{library_name}/log/{library_name}.log"}

    def _get_input_files_run(self, wildcards):
        inputs = reads_input_files(self.parent, wildcards.library_name)
        tool_config = self.config.get(self.name)
        for key in self.file_keys:
            if value := getattr(tool_config, key):
                inputs[key] = value
        return inputs

    @dictify
    def _get_output_files_run(self):
        # The tool writes one file per input file, named by params: directories, not files
        return (
            ("out", directory("work/{library_name}/out")),
            ("report", directory("work/{library_name}/report")),
            ("rejected", directory("work/{library_name}/rejected")),
        )

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
    wrappers = {"run": "bbduk/run"}

    @property
    def default_resource_usage(self) -> ResourceUsage:
        return ResourceUsage(threads=self.config.bbduk.num_threads, runtime="12h", mem="24000MB")

    file_keys = ("adapter_sequences", "barcodes")


class FastpStepPart(AdapterTrimmingStepPart):
    name = "fastp"
    wrappers = {"run": "fastp/run"}

    #: The HTML report of each FASTQ file or pair, below the path of the file
    report_outputs = {
        "run": {
            "report": ReportOutput("fastp", caption="fastp.rst", patterns=("{reads}/report.html",))
        }
    }

    @property
    def default_resource_usage(self) -> ResourceUsage:
        return ResourceUsage(threads=self.config.fastp.num_threads, runtime="12h", mem="24000MB")

    file_keys = ("filter_by_index1", "filter_by_index2")


class AdapterTrimmingWorkflow(BaseStep):
    name = "adapter_trimming"
    produces = [DataSignature(DataType.RAW, frozenset({"trimmed"}))]
    config_model_class = AdapterTrimmingConfigModel
    tool_parts = {Tool.BBDUK: BbdukStepPart, Tool.FASTP: FastpStepPart}

    @classmethod
    def get_output_paths(cls, config, signature=None, **kwargs) -> dict[str, str]:
        """Return local output paths for trimmed/raw FASTQ consumption."""
        return {"fastq_dir": "work"}
