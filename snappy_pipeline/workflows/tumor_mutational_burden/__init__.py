from biomedsheets.shortcuts import CancerCaseSheet, CancerCaseSheetOptions

from snappy_pipeline.utils import dictify
from snappy_pipeline.workflows.abstract import (
    BaseStep,
    BaseStepPart,
    ResourceUsage,
)
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType
from snappy_pipeline.workflows.variant_calling.model import ExpectedSomaticVariants

from .model import TumorMutationalBurden as TumorMutationalBurdenConfigModel


class TumorMutationalBurdenCalculationStepPart(BaseStepPart):
    """Calculation tumor mutational burden for each sample"""

    name = "tmb_gathering"

    actions = ("run",)

    entity_wildcard = "tumor_library"

    published_keys = ("json", "json_md5")

    default_resource_usage = ResourceUsage(threads=2, runtime="1h", mem="4096MB")

    @dictify
    def _get_input_files_run(self, wildcards):
        variants: ExpectedSomaticVariants = self.parent.get_upstream_paths(
            "variants", library_name=wildcards.tumor_library
        )
        yield "vcf", variants.vcf
        yield "vcf_tbi", variants.vcf_tbi
        yield "target_regions", self.config.target_regions

    def _get_output_files_run(self):
        tpl = "work/{tumor_library}/out/tmb.{tumor_library}.json"
        return {"json": tpl, "json_md5": tpl + ".md5"}

    def _get_log_file_run(self):
        return {"log": "work/{tumor_library}/log/tmb.{tumor_library}.log"}

    def _get_params_run(self, wildcards):
        return {
            "missense_re": self.config.missense_regex,
            "has_annotation": self.config.has_annotation,
        }


class TumorMutationalBurdenCalculationWorkflow(BaseStep):
    """Perform TMB calculation"""

    name = "tumor_mutational_burden"
    produces = [DataSignature(DataType.TABULAR, frozenset({"tmb"}))]
    config_model_class = TumorMutationalBurdenConfigModel
    sheet_shortcut_class = CancerCaseSheet
    sheet_shortcut_kwargs = {
        "options": CancerCaseSheetOptions(allow_missing_normal=True, allow_missing_tumor=True)
    }

    @classmethod
    def get_output_paths(cls, config, signature=None, **kwargs) -> dict[str, str]:
        """Return local TMB output paths for downstream consumers."""
        lib = kwargs.get("library_name", "{library_name}")
        return {"json": f"work/{lib}/out/tmb.{lib}.json"}

    def __init__(self, workflow, project, task_name):
        super().__init__(workflow, project, task_name)

        # Register sub step classes so the sub steps are available
        self.register_sub_step_classes((TumorMutationalBurdenCalculationStepPart,))
