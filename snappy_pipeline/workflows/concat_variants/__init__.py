"""Concatenate variant calls

Concatenates the VCFs of several tasks into one VCF per entity, for example the small variants and the structural variants of each pedigree after FDR control, so that one annotation, one filtration and one export follow.
``depends_on.variants`` names the tasks.
They must write the same entities, with the same ``group_by`` and ``library_selection``, and the task takes their grouping.
``bcftools concat`` also needs the same samples in the same order in every input, which the calls of one pedigree or one library satisfy.
The tool is ``bcftools concat --allow-overlaps``, which keeps every record of every input; ``extra_arguments`` such as ``--rm-dups exact`` drop records that two inputs share.
The output carries the tags that all inputs have, such as ``germline``, ``varlociraptor`` or ``filtered``, and the variant kinds of any of them, such as ``snv``, ``indel`` and ``sv``.

Outputs
-------

Below ``tasks/<task>/output/``:

- ``{library_name}/out/{library_name}.vcf.gz`` with ``.tbi``, ``.md5`` and ``.tbi.md5`` files.

- with ``calls_report``: ``{library_name}/report/calls/{library_name}.calls.tsv``, the calls of the VCF as a table (one row per call), and the same table as a datavzrd site below ``report/calls/html/``, which ``snappy report`` shows.

Downstream tasks read ``vcf`` and ``vcf_tbi`` through ``depends_on``.
"""

from typing import Any

from snakemake.iocontainers import Wildcards

from snappy_pipeline.base import ConfigError
from snappy_pipeline.utils import dictify
from snappy_pipeline.workflows.abstract import (
    BaseStep,
    BaseStepPart,
    ResourceUsage,
)
from snappy_pipeline.workflows.common import vcf_outputs
from snappy_pipeline.workflows.common.calls_report import CallsReportStepPart
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType

from .model import ConcatVariants as ConcatVariantsConfigModel

_OUT_PREFIX = "work/{library_name}/out/{library_name}"

#: Tags that name a kind of variant; the output carries those of any input
VARIANT_KINDS = frozenset({"snv", "indel", "mnv", "sv", "cnv", "mei", "repeats", "fusion"})


class BcftoolsConcatStepPart(BaseStepPart):
    """``bcftools concat`` over the VCFs of the upstream tasks, per entity"""

    name = "bcftools"
    actions = ("run",)
    wrappers = {"run": "bcftools/concat_variants"}
    result_keys = {"run": tuple(vcf_outputs(_OUT_PREFIX))}
    default_resource_usage = ResourceUsage(threads=1, mem="2G", runtime="2h")

    @dictify
    def _get_input_files_run(self, wildcards: Wildcards):
        vcf, tbi = [], []
        for task in self.parent.upstream_task_names("variants"):
            paths = self.parent.get_upstream_paths(
                "variants", upstream_task=task, library_name=wildcards.library_name
            )
            vcf.append(paths.vcf)
            tbi.append(paths.vcf_tbi)
        yield "vcf", vcf
        yield "vcf_tbi", tbi

    def _get_output_files_run(self) -> dict[str, str]:
        return vcf_outputs(_OUT_PREFIX)

    def _get_log_file_run(self) -> dict[str, str]:
        return {"log": "work/{library_name}/log/{library_name}.log"}

    def _get_params_run(self, wildcards: Wildcards) -> dict[str, Any]:
        return {}  # the base class adds the extra_arguments of the tool section


class ConcatVariantsWorkflow(BaseStep):
    name = "concat_variants"
    produces = [DataSignature(DataType.VARIANTS)]
    config_model_class = ConcatVariantsConfigModel

    @classmethod
    def task_produces(cls, config, upstream):
        """Concatenated calls: the tags that all inputs have, plus the variant kinds of any input"""
        signatures = upstream["variants"]
        common = frozenset.intersection(*(s.tags for s in signatures))
        kinds = frozenset(t for s in signatures for t in s.tags if t in VARIANT_KINDS)
        return (DataSignature(DataType.VARIANTS, common | kinds),)

    @classmethod
    def get_output_paths(cls, config, signature=None, **kwargs) -> dict[str, str]:
        lib = kwargs.get("library_name", "{library_name}")
        paths = vcf_outputs(_OUT_PREFIX.format(library_name=lib))
        return {key: paths[key] for key in ("vcf", "vcf_tbi")}

    def __init__(self, workflow, project, task_name):
        super().__init__(workflow, project, task_name)
        self.register_sub_step_classes((BcftoolsConcatStepPart,))
        if self.config.calls_report:
            self.register_sub_step_classes((CallsReportStepPart,))

    def check_config(self):
        """The inputs must write the same entities, so every input task has to declare the
        grouping of this task."""
        group_by = self.effective_group_by
        for task in self.upstream_task_names("variants"):
            declared = self.project.group_by[task]
            if declared != group_by:
                raise ConfigError(
                    f"the task concatenates with group_by: {group_by}, but its input task "
                    f"{task!r} writes with group_by: {declared}"
                )
