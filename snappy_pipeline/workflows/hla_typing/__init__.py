# -*- coding: utf-8 -*-
"""HLA typing from sequencing reads

The step determines the HLA type of each library from its reads.
The ``tool`` field selects the typer:

- ``optitype``: types DNA and RNA libraries from the FASTQ files.
  It reads ``depends_on.reads``, which names ``data_sets``, an ``external_data`` task, or an ``adapter_trimming`` task (:ref:`step_adapter_trimming`).
- ``arcashla``: types RNA-seq libraries from the BAM files of ``depends_on.alignments``, for example an :ref:`step_ngs_mapping` task with ``tool: star``.
  A task takes all RNA libraries unless it sets ``library_selection``.
  A library is single-end if no read of its BAM file is paired.
  The options ``avg`` and ``std`` give the mean and standard deviation of the fragment length, for single-end libraries only.

The result is one set of files per library.

Outputs
-------

For each library ``{library_name}``, the task links the following result files below ``tasks/<task>/output/``.
Each file has an ``.md5`` file.

- ``{library_name}/out/{library_name}.txt``
- ``{library_name}/out/{library_name}.json``

With ``optitype``, the result table and the coverage plot are result files, too:

- ``{library_name}/out/{library_name}.result.tsv``
- ``{library_name}/out/{library_name}.coverage_plot.pdf``

Downstream tasks read the keys ``txt`` and ``calls_json`` (the ``.json`` file).
"""

import re

from snappy_pipeline.utils import dictify
from snappy_pipeline.workflows.abstract import (
    RuleOptions,
    BaseStep,
    BaseStepPart,
    ReportOutput,
    ResourceUsage,
)
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType
from snappy_pipeline.workflows.common.reads import reads_input_files, reads_params

from .model import HlaTyping as HlaTypingConfigModel
from .model import Tool

#: Output keys and extensions of the main payload
_OUTPUTS = {"txt": ".txt", "txt_md5": ".txt.md5", "json": ".json", "json_md5": ".json.md5"}

#: Output keys and extensions that only ``optitype`` writes
_OPTITYPE_OUTPUTS = {
    "tsv": ".result.tsv",
    "tsv_md5": ".result.tsv.md5",
    "cov_pdf": ".coverage_plot.pdf",
    "cov_pdf_md5": ".coverage_plot.pdf.md5",
}

_OUT_PATTERN = "work/{library_name}/out/{library_name}"


class OptiTypeStepPart(BaseStepPart):
    """HLA Typing using OptiType"""

    #: Step name
    name = "optitype"

    #: Class available actions
    actions = ("run",)

    wrappers = {"run": "optitype"}
    result_keys = tuple(_OUTPUTS | _OPTITYPE_OUTPUTS)
    report_outputs = {"run": {"cov_pdf": ReportOutput("OptiType", caption="optitype.rst")}}

    @property
    def default_resource_usage(self) -> ResourceUsage:
        threads = self.config.optitype.num_mapping_threads
        return ResourceUsage(threads=threads, runtime="40h", mem="45000MB")

    #: The MHC classes that the tool types
    mhc_classes = ("class_i",)

    supported_extraction_types = ("dna", "rna")

    def _get_input_files_run(self, wildcards):
        """Return input files"""
        return reads_input_files(self.parent, wildcards.library_name)

    def _get_output_files_run(self):
        outputs = _OUTPUTS | _OPTITYPE_OUTPUTS
        return {key: _OUT_PATTERN + ext for key, ext in outputs.items()}

    def _get_log_file_run(self):
        return {"log": "work/{library_name}/log/{library_name}.log"}

    def _get_params_run(self, wildcards):
        """Return dict for the wrapper"""
        result = reads_params(self.parent, wildcards.library_name)
        result["seq_type"] = self._get_seq_type(wildcards)
        result["use_discordant"] = "true" if self.config.optitype.use_discordant else "false"
        result["num_mapping_threads"] = self.config.optitype.num_mapping_threads
        result["solver"] = self.config.optitype.solver
        result["unpaired_weight"] = self.config.optitype.unpaired_weight
        result["max_reads"] = self.config.optitype.max_reads
        result["yara_error_rate"] = self.config.optitype.yara_mapper.error_rate
        result["yara_strata_rate"] = self.config.optitype.yara_mapper.strata_rate
        result["yara_sensitivity"] = self.config.optitype.yara_mapper.sensitivity
        return result

    def _get_seq_type(self, wildcards):
        """Return sequence type for the library name in wildcards"""
        return self.parent.library_row(wildcards.library_name)["extraction_type"]


class ArcasHlaStepPart(BaseStepPart):
    """HLA Typing using arcasHLA"""

    #: Step name
    name = "arcashla"

    #: Class available actions
    actions = ("prepare_reference", "run")

    wrappers = {"prepare_reference": "arcashla/prepare_reference", "run": "arcashla/run"}
    rule_options = {"prepare_reference": RuleOptions(retries=2)}
    extra_arguments_action = "run"
    result_keys = {"run": tuple(_OUTPUTS)}
    default_resource_usage = ResourceUsage(threads=4, runtime="60h", mem="64000MB")

    #: Downloading and indexing the reference
    resource_usage = {"prepare_reference": ResourceUsage(threads=1, runtime="4h", mem="8GB")}

    #: The MHC classes that the tool types
    mhc_classes = ("class_i", "class_ii")

    supported_extraction_types = ("rna",)

    PAIRED_PATTERN: re.Pattern = re.compile(
        r"^(?P<passPairs>[0-9]+) \+ (?P<failPairs>[0-9]+) paired in sequencing$"
    )

    def _get_output_files_prepare_reference(self):
        return {"done": "work/prepare_reference/out/.done"}

    def _get_log_file_prepare_reference(self):
        return {"log": "work/prepare_reference/log/prepare_reference.log"}

    def _get_params_prepare_reference(self, wildcards):
        return {"reference_version": self.config.arcashla.reference_version}

    @dictify
    def _get_input_files_run(self, wildcards):
        """Return input files"""
        yield "ref_done", "work/prepare_reference/out/.done"
        alignments = self.parent.get_upstream_paths(
            "alignments", library_name=wildcards.library_name
        )
        yield "bam", alignments.bam

    def _get_output_files_run(self):
        return {key: _OUT_PATTERN + ext for key, ext in _OUTPUTS.items()}

    def _get_log_file_run(self):
        return {"log": "work/{library_name}/log/{library_name}.log"}

    def _get_params_run(self, wildcards):
        """The options of ``arcasHLA genotype``; ``avg`` and ``std`` are for single-end reads"""
        return self.config.arcashla.model_dump(
            mode="json", exclude={"extra_arguments", "reference_version"}, exclude_none=True
        )


class HlaTypingWorkflow(BaseStep):
    """Perform HLA Typing"""

    #: Step name
    name = "hla_typing"

    produces = [DataSignature(DataType.TABULAR, frozenset({"hla"}))]

    #: config_model_class
    config_model_class = HlaTypingConfigModel

    tool_parts = {Tool.optitype: OptiTypeStepPart, Tool.arcashla: ArcasHlaStepPart}

    @classmethod
    def task_library_requirement(cls, config):
        """The libraries whose extraction type the tool supports."""
        return f"extraction_type in {list(cls.tool_parts[config.tool].supported_extraction_types)}"

    @classmethod
    def task_library_default(cls, config):
        """arcashla takes all RNA libraries, which are all it supports."""
        return "" if config.tool == "arcashla" else None

    @classmethod
    def get_output_paths(cls, config, signature=None, **kwargs) -> dict[str, str]:
        """Return local HLA typing output paths for downstream consumers."""
        lib = kwargs.get("library_name", "{library_name}")
        return {"txt": f"work/{lib}/out/{lib}.txt", "calls_json": f"work/{lib}/out/{lib}.json"}

    @classmethod
    def mhc_classes(cls, config) -> tuple[str, ...]:
        """The MHC classes that a task of this step types."""
        return cls.tool_parts[config.tool].mhc_classes
