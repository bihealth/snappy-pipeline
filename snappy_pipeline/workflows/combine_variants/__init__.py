"""Combine somatic and germline variants

Merges the somatic VCF of a tumor library with the germline VCF of its matched normal library.
The task needs both ``depends_on.somatic_variants`` and ``depends_on.germline_variants``.
It produces one VCF per tumor library, which needs a matched normal library.
The merge runs GATK 3 ``CombineVariants``.

With ``rename_combined`` set to ``tumor`` or ``germline``, the sample in the combined VCF is named after the tumor or the matched normal library.

Outputs
-------

Below ``tasks/<task>/output/``:

- ``{library_name}/out/{library_name}.vcf.gz`` with ``.tbi``, ``.md5`` and ``.tbi.md5`` files.

Downstream tasks read ``vcf`` and ``vcf_tbi`` through ``depends_on``.
"""

from typing import Any

from biomedsheets.shortcuts import CancerCaseSheet
from snakemake.io import expand
from snakemake.iocontainers import Wildcards

from snappy_pipeline.base import MissingConfiguration
from snappy_pipeline.models import RelationshipDefinition
from snappy_pipeline.utils import dictify
from snappy_pipeline.workflows.abstract import (
    BaseStep,
    BaseStepPart,
    ResourceUsage,
)
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType

from .model import CombineVariants as CombineVariantsConfigModel
from .model import RenameCombine

__author__ = "Eric Blanc <eric.blanc@bih-charite.de>"


_OUT_PREFIX = "work/{tumor_library}/out/{tumor_library}"
_LOG_PREFIX = "work/{tumor_library}/log/{tumor_library}"

#: Output keys and extensions of the combined VCF
_OUTPUTS = {
    "vcf": ".vcf.gz",
    "vcf_tbi": ".vcf.gz.tbi",
    "vcf_md5": ".vcf.gz.md5",
    "vcf_tbi_md5": ".vcf.gz.tbi.md5",
}


class CombineVariantsStepPart(BaseStepPart):
    name = "combine"
    actions = ("run",)
    default_resource_usage = ResourceUsage(threads=1, mem="4G", runtime="4h")

    @dictify
    def _get_input_files_run(self, wildcards: Wildcards):
        yield from self.parent.reference_inputs().items()

        somatic = self.parent.get_upstream_paths(
            "somatic_variants", library_name=wildcards.tumor_library
        )
        yield "somatic_vcf", somatic.vcf

        normal_lib = self.get_normal_lib_name(wildcards.tumor_library)
        germline = self.parent.get_upstream_paths("germline_variants", library_name=normal_lib)
        yield "germline_vcf", germline.vcf

    def get_normal_lib_name(self, tumor_library: str) -> str:
        """Return the matched normal library of the tumor library; fail if there is none"""
        df = self.parent.build_library_dataframe()
        tumor_df = df[df["library_name"] == tumor_library]
        normal_lib = None if tumor_df.empty else tumor_df.iloc[0].get("matched_normal_lib")
        if not normal_lib:
            raise ValueError(
                f"Tumor library {tumor_library!r} has no matched normal library, "
                "which combine_variants needs for the germline variants"
            )
        return normal_lib

    def _get_output_files_run(self) -> dict[str, str]:
        return {key: f"{_OUT_PREFIX}{ext}" for key, ext in _OUTPUTS.items()}

    def _get_params_run(self, wildcards: Wildcards) -> dict[str, Any]:
        if name_origin := self.config.rename_combined:
            if name_origin == RenameCombine.TUMOR:
                sample_name = wildcards.tumor_library
            elif name_origin == RenameCombine.GERMLINE:
                sample_name = self.get_normal_lib_name(wildcards.tumor_library)
            else:
                raise MissingConfiguration(f"Unimplemented sample name type {name_origin}")
            return {"tumor_library": wildcards.tumor_library, "sample_name": sample_name}
        else:
            return {"tumor_library": wildcards.tumor_library}

    def get_result_files(self) -> list[str]:
        df = self.parent.build_library_dataframe()
        result = []
        for tumor_library in df[df["role"] == "tumor"]["library_name"].unique():
            self.get_normal_lib_name(tumor_library)
            for path in self.get_output_files("run").values():
                result += expand(path.replace("work/", "output/", 1), tumor_library=[tumor_library])
        return result

    def _get_log_file_run(self) -> dict[str, str]:
        return {"log": _LOG_PREFIX + ".log"}


class CombineVariantsWorkflow(BaseStep):
    name = "combine_variants"
    sheet_shortcut_class = CancerCaseSheet

    produces = [DataSignature(DataType.VARIANTS, frozenset({"combined"}))]

    config_model_class = CombineVariantsConfigModel

    default_relationships = {
        "matched_normal_lib": RelationshipDefinition(
            via="donor_name",
            target="role == 'normal' and extraction_type == 'dna'",
        )
    }

    @classmethod
    def get_output_paths(cls, config, signature=None, **kwargs) -> dict[str, str]:
        lib = kwargs.get("library_name", "{library_name}")
        return {
            "vcf": f"work/{lib}/out/{lib}.vcf.gz",
            "vcf_tbi": f"work/{lib}/out/{lib}.vcf.gz.tbi",
        }

    def __init__(self, workflow, project, task_name):

        super().__init__(workflow, project, task_name)
        self.register_sub_step_classes((CombineVariantsStepPart,))
