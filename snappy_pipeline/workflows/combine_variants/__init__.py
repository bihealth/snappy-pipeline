"""Combine somatic and germline variants

Merges the somatic VCF of a tumor library with the germline VCF of its matched normal library.
The task needs both ``depends_on.somatic_variants`` and ``depends_on.germline_variants``.
It writes one VCF per tumor DNA library.
A tumor DNA library without a matched normal library is an error.
``library_selection`` can only narrow the set of tumor DNA libraries.
The merge runs GATK 3 ``CombineVariants``.
The somatic sites carry the ``SOMATIC`` flag in the combined VCF.
The output carries the tags that both inputs have.
So if both inputs are annotated, the combined calls are annotated, too.
Neoepitope prediction, for example, needs calls annotated by VEP.

``rename_combined`` names the sample in the combined VCF after the tumor library (``tumor``) or the matched normal library (``germline``).

Outputs
-------

Below ``tasks/<task>/output/``:

- ``{library_name}/out/{library_name}.vcf.gz`` with ``.tbi``, ``.md5`` and ``.tbi.md5`` files.

- with ``calls_report``: ``{library_name}/report/calls/{library_name}.calls.tsv``, the calls of the VCF as a table with one row per call.
  ``snappy report`` shows the same table as a datavzrd site below ``report/calls/html/``.

Downstream tasks read ``vcf`` and ``vcf_tbi`` through ``depends_on``.
"""

from typing import Any

from snakemake.io import expand
from snakemake.iocontainers import Wildcards

from snappy_pipeline.base import ConfigError
from snappy_pipeline.utils import dictify
from snappy_pipeline.workflows.abstract import (
    MATCHED_NORMAL_RELATIONSHIP,
    BaseStep,
    BaseStepPart,
    ResourceUsage,
)
from snappy_pipeline.workflows.common import vcf_outputs
from snappy_pipeline.workflows.common.calls_report import CallsReportStepPart
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType

from .model import CombineVariants as CombineVariantsConfigModel
from .model import RenameCombine

__author__ = "Eric Blanc <eric.blanc@bih-charite.de>"


_OUT_PREFIX = "work/{tumor_library}/out/{tumor_library}"
_LOG_PREFIX = "work/{tumor_library}/log/{tumor_library}"


class CombineVariantsStepPart(BaseStepPart):
    name = "combine"
    actions = ("run",)
    wrappers = {"run": "combine_variants/run"}
    default_resource_usage = ResourceUsage(threads=1, mem="4G", runtime="4h")

    def tool_section(self):
        """The options of GATK's CombineVariants are on the step model"""
        return self.config

    def get_result_files(self) -> list[str]:
        result = []
        for tumor_library in self.parent.libraries["library_name"].unique():
            self.get_normal_lib_name(tumor_library)
            for path in self.get_output_files("run").values():
                result += expand(path, tumor_library=[tumor_library])
        return result

    @dictify
    def _get_input_files_run(self, wildcards: Wildcards):
        yield from self.parent.reference_inputs().items()

        somatic = self.parent.get_upstream_paths(
            "somatic_variants", library_name=wildcards.tumor_library
        )
        yield "somatic_vcf", somatic.vcf

        normal_lib = self.get_normal_lib_name(wildcards.tumor_library)
        entity = self.parent.upstream_entity("germline_variants", normal_lib)
        germline = self.parent.get_upstream_paths("germline_variants", library_name=entity)
        yield "germline_vcf", germline.vcf

    def _get_output_files_run(self) -> dict[str, str]:
        return vcf_outputs(_OUT_PREFIX)

    def _get_log_file_run(self) -> dict[str, str]:
        return {"log": _LOG_PREFIX + ".log"}

    def _get_params_run(self, wildcards: Wildcards) -> dict[str, Any]:
        if name_origin := self.config.rename_combined:
            if name_origin == RenameCombine.TUMOR:
                sample_name = wildcards.tumor_library
            elif name_origin == RenameCombine.GERMLINE:
                sample_name = self.get_normal_lib_name(wildcards.tumor_library)
            else:
                raise NotImplementedError(f"Unimplemented sample name type {name_origin}")
            return {"tumor_library": wildcards.tumor_library, "sample_name": sample_name}
        else:
            return {"tumor_library": wildcards.tumor_library}

    def get_normal_lib_name(self, tumor_library: str) -> str:
        """Return the matched normal library of the tumor library; fail if there is none"""
        normal_lib = self.parent.matched_normal(tumor_library)
        if not normal_lib:
            raise ConfigError(
                f"Tumor library {tumor_library!r} has no matched normal library, "
                "which combine_variants needs for the germline variants"
            )
        return normal_lib


class CombineVariantsWorkflow(BaseStep):
    name = "combine_variants"
    produces = [DataSignature(DataType.VARIANTS, frozenset({"combined"}))]
    config_model_class = CombineVariantsConfigModel

    #: The tumor DNA libraries, whose somatic and germline calls are combined
    library_requirement = "role == 'tumor' and extraction_type == 'dna'"

    default_relationships = MATCHED_NORMAL_RELATIONSHIP

    @classmethod
    def task_produces(cls, config, upstream):
        """Combined calls, with the tags that both inputs have, such as ``annotated`` and the
        annotator."""
        common = upstream["somatic_variants"].tags & upstream["germline_variants"].tags
        return (cls.produces[0].with_tags(*common),)

    @classmethod
    def get_output_paths(cls, config, signature=None, **kwargs) -> dict[str, str]:
        lib = kwargs.get("library_name", "{library_name}")
        paths = vcf_outputs(_OUT_PREFIX.format(tumor_library=lib))
        return {key: paths[key] for key in ("vcf", "vcf_tbi")}

    @classmethod
    def default_group_by(cls, config):
        """One combined VCF per tumor library, whatever grouping the germline calls have."""
        return "library"

    def __init__(self, workflow, project, task_name):

        super().__init__(workflow, project, task_name)
        self.register_sub_step_classes((CombineVariantsStepPart,))
        if self.config.calls_report:
            self.register_sub_step_classes((CallsReportStepPart,))
