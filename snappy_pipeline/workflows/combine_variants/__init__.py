from typing import Any

from biomedsheets.shortcuts import CancerCaseSheet
from snakemake.iocontainers import Wildcards

from snappy_pipeline.base import MissingConfiguration
from snappy_pipeline.models import RelationshipDefinition
from snappy_pipeline.utils import dictify, listify
from snappy_pipeline.workflows.abstract import (
    BaseStep,
    BaseStepPart,
    LinkOutStepPart,
    ResourceUsage,
)
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType
from snappy_pipeline.workflows.variant_annotation import VariantAnnotationWorkflow
from snappy_pipeline.workflows.variant_calling import VariantCallingWorkflow
from snappy_pipeline.workflows.variant_filtration import VariantFiltrationWorkflow

from .model import CombineVariants as CombineVariantsConfigModel
from .model import RenameCombine

__author__ = "Eric Blanc <eric.blanc@bih-charite.de>"

DEFAULT_CONFIG = CombineVariantsConfigModel.default_config_yaml_string()

_OUT_PREFIX = "work/{tumor_library}/out/{tumor_library}"
_LOG_PREFIX = "work/{tumor_library}/log/{tumor_library}"


class CombineVariantsStepPart(BaseStepPart):
    name = "combine"
    actions = ("run",)
    default_resource_usage = ResourceUsage(threads=1, mem="4G", runtime="4h")

    def get_input_files(self, action: str):
        self._validate_action(action)
        return getattr(self, f"_get_input_files_{action}")

    @dictify
    def _get_input_files_run(self, wildcards: Wildcards):
        yield "reference", self.w_config.static_data_config.reference.path

        somatic = self.parent.get_upstream_paths(
            "somatic_variant", library_name=wildcards.tumor_library
        )
        yield "somatic_vcf", somatic.vcf

        df = self.parent.build_library_dataframe()
        tumor_df = df[df["library_name"] == wildcards.tumor_library]
        normal_lib = None
        if not tumor_df.empty:
            normal_lib = tumor_df.iloc[0].get("matched_normal_lib") or None
        if normal_lib:
            germline = self.parent.get_upstream_paths("germline_variant", library_name=normal_lib)
            yield "germline_vcf", germline.vcf

    def get_output_files(self, action: str) -> dict[str, Any]:
        match action:
            case "run":
                return {
                    k: _OUT_PREFIX + ext
                    for k, ext in (
                        ("vcf", ".vcf.gz"),
                        ("vcf_tbi", ".vcf.gz.tbi"),
                        ("vcf_md5", ".vcf.gz.md5"),
                        ("vcf_tbi_md5", ".vcf.gz.tbi.md5"),
                    )
                }
            case _:
                raise MissingConfiguration(f"Unimplemented action {action}")

    def get_args(self, action: str):
        self._validate_action(action)
        return getattr(self, f"_get_args_{action}")

    def _get_args_run(self, wildcards: Wildcards) -> dict[str, Any]:
        if name_origin := self.config.rename_combined:
            df = self.parent.build_library_dataframe()
            tumor_df = df[df["library_name"] == wildcards.tumor_library]
            if name_origin == RenameCombine.TUMOR:
                sample_name = wildcards.tumor_library
            elif name_origin == RenameCombine.GERMLINE:
                if not tumor_df.empty:
                    sample_name = tumor_df.iloc[0].get("matched_normal_lib") or ""
                else:
                    sample_name = ""
            else:
                raise MissingConfiguration(f"Unimplemented sample name type {name_origin}")
            return {"tumor_library": wildcards.tumor_library, "sample_name": sample_name}
        else:
            return {"tumor_library": wildcards.tumor_library}

    @dictify
    def get_log_file(self, action: str):
        match action:
            case "run":
                for k, ext in (
                    ("log", ".log"),
                    ("conda_list", ".conda_list.txt"),
                    ("conda_info", ".conda_info.txt"),
                ):
                    yield k, _LOG_PREFIX + ext
                    yield k + "_md5", _LOG_PREFIX + ext + ".md5"
            case _:
                raise MissingConfiguration(f"Unimplemented action {action}")


class CombineVariantsWorkflow(BaseStep):
    name = "combine_variants"
    sheet_shortcut_class = CancerCaseSheet

    consumes = {DataSignature(DataType.VARIANTS): True}
    produces = [DataSignature(DataType.VARIANTS, frozenset({"combined"}))]

    config_model_class = CombineVariantsConfigModel

    default_relationships = {
        "matched_normal_lib": RelationshipDefinition(
            via="donor_name",
            target="role == 'normal' and extraction_type == 'dna'",
        )
    }

    @classmethod
    def default_config_yaml(cls):
        return DEFAULT_CONFIG

    @classmethod
    def get_output_paths(cls, signature=None, **kwargs) -> dict[str, str]:
        cls.require_signature(signature)
        lib = kwargs.get("library_name", "{library_name}")
        return {
            "vcf": f"output/{lib}/out/{lib}.vcf.gz",
            "vcf_tbi": f"output/{lib}/out/{lib}.vcf.gz.tbi",
        }

    def __init__(self, workflow, config, config_lookup_paths, config_paths, workdir, **kwargs):
        previous_steps = [
            VariantCallingWorkflow,
            VariantAnnotationWorkflow,
            VariantFiltrationWorkflow,
        ]

        super().__init__(
            workflow,
            config,
            config_lookup_paths,
            config_paths,
            workdir,
            previous_steps=previous_steps,
            **kwargs,
        )
        self.register_sub_step_classes((CombineVariantsStepPart, LinkOutStepPart))

    @listify
    def get_result_files(self):
        df = self.build_library_dataframe()
        tumor_libs = df[df["role"] == "tumor"]["library_name"].unique()
        for t in tumor_libs:
            for ext in ("", ".tbi", ".md5", ".tbi.md5"):
                yield f"output/{t}/out/{t}.vcf.gz{ext}"
            for ext in ("log", "conda_list.txt", "conda_info.txt"):
                for hash_suffix in ("", ".md5"):
                    yield f"output/{t}/log/{t}.{ext}{hash_suffix}"
