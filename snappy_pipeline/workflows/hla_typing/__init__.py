# -*- coding: utf-8 -*-
"""HLA typing from sequencing reads

The step determines the HLA type of each library from its reads.
The ``tool`` field selects the typer:

- ``optitype``: types DNA and RNA libraries from the FASTQ files. It reads ``depends_on.reads``, which names ``data_sets``, an ``external_data`` task, or an ``adapter_trimming`` task (:ref:`step_adapter_trimming`).
- ``arcashla``: types RNA-seq libraries from the BAM files of ``depends_on.alignments``, for example an :ref:`step_ngs_mapping` task with ``tool: star``.

The result is one set of files per library.

Outputs
-------

The task publishes the following files below ``tasks/<task>/output/``, for each library with name ``{library_name}``.
Each file has an ``.md5`` file.

- ``{library_name}/out/{library_name}.txt``
- ``{library_name}/out/{library_name}.json``

Downstream tasks read the keys ``txt`` and ``calls_json`` (the ``.json`` file).
With ``optitype``, the result table and the coverage plot stay in ``work/{library_name}/out/``.
"""

import re
from collections import OrderedDict

from biomedsheets.shortcuts import GenericSampleSheet
from snakemake.io import expand

from snappy_pipeline.utils import dictify, listify
from snappy_pipeline.workflows.abstract import (
    BaseStep,
    BaseStepPart,
    ResourceUsage,
)
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType
from snappy_pipeline.workflows.common.reads import reads_input_files, reads_params

from .model import HlaTyping as HlaTypingConfigModel

#: Output keys and extensions of the main payload
_OUTPUTS = {"txt": ".txt", "txt_md5": ".txt.md5", "json": ".json", "json_md5": ".json.md5"}

_OUT_PATTERN = "work/{library_name}/out/{library_name}"

#: HLA typing tools
HLA_TYPERS = ("optitype", "arcashla")

#: Default configuration for the hla_typing schema


class OptiTypeStepPart(BaseStepPart):
    """HLA Typing using OptiType"""

    #: Step name
    name = "optitype"

    supported_extraction_types = ("dna", "rna")

    #: Class available actions
    actions = ("run",)

    def _get_input_files_run(self, wildcards):
        """Return input files"""
        return reads_input_files(self.parent, wildcards.library_name)

    def _get_output_files_run(self):
        if self.name != str(self.config.tool):
            return {}
        # additional optitype output files
        outputs = _OUTPUTS | {"tsv": ".result.tsv", "cov_pdf": ".coverage_plot.pdf"}
        return {key: _OUT_PATTERN + ext for key, ext in outputs.items()}

    def _get_log_file_run(self):
        if self.name != str(self.config.tool):
            return {}
        return {"log": "work/{library_name}/log/{library_name}.log"}

    def _get_params_run(self, wildcards):
        """Return dict for the wrapper"""
        result = reads_params(self.parent, wildcards.library_name)
        result["seq_type"] = self._get_seq_type(wildcards)
        result["use_discordant"] = "true" if self.config.optitype.use_discordant else "false"
        result["num_mapping_threads"] = self.config.optitype.num_mapping_threads
        result["max_reads"] = self.config.optitype.max_reads
        result["yara_error_rate"] = self.config.optitype.yara_mapper.error_rate
        result["yara_strata_rate"] = self.config.optitype.yara_mapper.strata_rate
        result["yara_sensitivity"] = self.config.optitype.yara_mapper.sensitivity
        return result

    def _get_seq_type(self, wildcards):
        """Return sequence type for the library name in wildcards"""
        library = self.parent.ngs_library_name_to_ngs_library[wildcards.library_name]
        return library.test_sample.extra_infos.get("extractionType", "DNA").lower()

    default_resource_usage = ResourceUsage(threads=6, runtime="40h", mem="45000MB")


class ArcasHlaStepPart(BaseStepPart):
    """HLA Typing using arcasHLA"""

    #: Step name
    name = "arcashla"

    supported_extraction_types = ("rna",)

    #: Class available actions
    actions = ("run",)

    PAIRED_PATTERN: re.Pattern = re.compile(
        r"^(?P<passPairs>[0-9]+) \+ (?P<failPairs>[0-9]+) paired in sequencing$"
    )

    def __init__(self, parent):
        super().__init__(parent)
        self.mapper = self.config.arcashla.mapper

    @dictify
    def _get_input_files_run(self, wildcards):
        """Return input files"""
        yield "ref_done", "work/prepare_reference/out/.done"
        alignments = self.parent.get_upstream_paths(
            "alignments", library_name=wildcards.library_name
        )
        yield "bam", alignments.bam

    def _get_output_files_run(self):
        if self.name != str(self.config.tool):
            return {}
        return {key: _OUT_PATTERN + ext for key, ext in _OUTPUTS.items()}

    def _get_log_file_run(self):
        return {"log": "work/{library_name}/log/{library_name}.log"}

    default_resource_usage = ResourceUsage(threads=4, runtime="60h", mem="15000MB")


class HlaTypingWorkflow(BaseStep):
    """Perform HLA Typing"""

    #: Step name
    name = "hla_typing"
    produces = [DataSignature(DataType.TABULAR, frozenset({"hla"}))]

    #: Default biomed sheet class
    sheet_shortcut_class = GenericSampleSheet

    #: config_model_class
    config_model_class = HlaTypingConfigModel

    @classmethod
    def get_output_paths(cls, config, signature=None, **kwargs) -> dict[str, str]:
        """Return local HLA typing output paths for downstream consumers."""
        lib = kwargs.get("library_name", "{library_name}")
        return {"txt": f"work/{lib}/out/{lib}.txt", "calls_json": f"work/{lib}/out/{lib}.json"}

    def __init__(self, workflow, project, task_name):
        super().__init__(workflow, project, task_name)
        match self.config.tool:
            case "optitype":
                selected = OptiTypeStepPart
            case "arcashla":
                selected = ArcasHlaStepPart
            case _:
                raise NotImplementedError(f"Unknown tool: {self.config.tool}")
        self.register_sub_step_classes((selected,))
        #: Mapping from library name to library object
        self.ngs_library_name_to_ngs_library = OrderedDict()
        for sheet in self.shortcut_sheets:
            for ngs_library in sheet.all_ngs_libraries:
                self.ngs_library_name_to_ngs_library[ngs_library.name] = ngs_library

    @listify
    def get_result_files(self):
        """Return list of result files for the NGS mapping workflow

        We will process all NGS libraries of all test samples in all sample
        sheets.
        """
        yield from self._yield_result_files(
            "output/{library_name}/out/{library_name}{ext}", ext=_OUTPUTS.values()
        )

    def _yield_result_files(self, tpl, **kwargs):
        """Build output paths from path template and extension list"""
        tool = str(self.config.tool)
        supported = self.sub_steps[tool].supported_extraction_types
        df = self.build_library_dataframe()
        if df.empty:
            return
        for library_name in self.output_entities:
            row = df[df["library_name"] == library_name]
            if row.empty:
                continue
            extraction_type = row.iloc[0].get("extraction_type", "").lower()
            if extraction_type in supported:
                yield from expand(
                    tpl,
                    library_name=[library_name],
                    **kwargs,
                )
