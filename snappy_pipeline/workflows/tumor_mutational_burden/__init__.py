import os

from biomedsheets.shortcuts import CancerCaseSheet, CancerCaseSheetOptions
from snakemake.io import expand

from snappy_pipeline.utils import dictify, listify
from snappy_pipeline.workflows.abstract import (
    BaseStep,
    BaseStepPart,
    LinkOutStepPart,
    ResourceUsage,
)
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType
from snappy_pipeline.workflows.variant_calling.model import ExpectedSomaticVariants

from .model import TumorMutationalBurden as TumorMutationalBurdenConfigModel

#: Extensions of files to create as main payload
EXT_VALUES = (".json", ".json.md5")

#: Names of the files to create for the extension
EXT_NAMES = ("json", "json_md5")

#: Default configuration for the tmb calculation step


class TumorMutationalBurdenCalculationStepPart(BaseStepPart):
    """Calculation tumor mutational burden for each sample"""

    name = "tmb_gathering"

    actions = ("run",)

    @dictify
    def get_input_files(self, action):
        self._validate_action(action)

        variants: ExpectedSomaticVariants = self.parent.get_upstream_paths(
            "somatic_variant", library_name="{tumor_library}"
        )
        yield "vcf", variants.vcf
        yield "vcf_tbi", variants.vcf_tbi

    @dictify
    def get_output_files(self, action):
        # Validate action
        self._validate_action(action)

        base_name = "tmb.{tumor_library}"

        tpl = os.path.join("work", "{tumor_library}", "out", base_name)

        key_ext = {"json": ".json"}
        for key, ext in key_ext.items():
            yield key, tpl + ext
            yield key + "_md5", tpl + ext + ".md5"

    @dictify
    def get_log_file(self, action):
        self._validate_action(action)

        base_name = "tmb.{tumor_library}"

        tpl = os.path.join("work", "{tumor_library}", "log", base_name)

        key_ext = (
            ("log", ".log"),
            ("conda_info", ".conda_info.txt"),
            ("conda_list", ".conda_list.txt"),
        )
        for key, ext in key_ext:
            yield key, tpl + ext
            yield key + "_md5", tpl + ext + ".md5"

    def get_resource_usage(self, action: str, **kwargs) -> ResourceUsage:
        self._validate_action(action)
        mem_mb = 4 * 1024  # 4GB
        return ResourceUsage(
            threads=2,
            runtime="1h",  # 1 hour
            mem=f"{mem_mb}MB",
        )

    def get_args(self, action):
        self._validate_action(action)
        return self._get_args_run

    def _get_args_run(self, _wildcards):
        return {
            "missense_re": self.config.missense_regex,
            "target_regions": self.config.target_regions,
            "has_annotation": self.config.has_annotation,
        }


class TumorMutationalBurdenCalculationWorkflow(BaseStep):
    """Perform TMB calculation"""

    name = "tumor_mutational_burden"
    consumes = {DataSignature(DataType.VARIANTS, frozenset({"somatic", ("snv", "indel")})): True}
    produces = [DataSignature(DataType.TABULAR, frozenset({"tmb"}))]
    config_model_class = TumorMutationalBurdenConfigModel
    sheet_shortcut_class = CancerCaseSheet
    sheet_shortcut_kwargs = {
        "options": CancerCaseSheetOptions(allow_missing_normal=True, allow_missing_tumor=True)
    }

    @classmethod
    def get_output_paths(cls, signature=None, **kwargs) -> dict[str, str]:
        """Return local TMB output paths for downstream consumers."""
        cls.require_signature(signature)
        lib = kwargs.get("library_name", "{library_name}")
        return {"json": f"output/{lib}/out/tmb.{lib}.json"}

    def __init__(
        self,
        workflow,
        config,
        config_lookup_paths,
        config_paths,
        workdir,
        task_name: str | None = None,
        **kwargs,
    ):
        super().__init__(
            workflow,
            config,
            config_lookup_paths,
            config_paths,
            workdir,
            previous_steps=(),
            task_name=task_name,
            **kwargs,
        )

        # Register sub step classes so the sub steps are available
        self.register_sub_step_classes((TumorMutationalBurdenCalculationStepPart, LinkOutStepPart))

    @listify
    def get_result_files(self):
        """Return list of result files for the TMB workflow."""
        log_exts = (
            ".log",
            ".log.md5",
            ".conda_info.txt",
            ".conda_info.txt.md5",
            ".conda_list.txt",
            ".conda_list.txt.md5",
        )
        for entity in self.output_entities:
            yield from expand(
                os.path.join("output", "{tumor_library}", "out", "tmb.{tumor_library}{ext}"),
                tumor_library=[entity],
                ext=EXT_VALUES,
            )
            yield from expand(
                os.path.join("output", "{tumor_library}", "log", "tmb.{tumor_library}{ext}"),
                tumor_library=[entity],
                ext=log_exts,
            )
