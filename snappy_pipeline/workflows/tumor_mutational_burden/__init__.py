import os
import sys
from collections import OrderedDict

from biomedsheets.shortcuts import CancerCaseSheet, CancerCaseSheetOptions, is_not_background
from snakemake.io import expand

from snappy_pipeline.utils import dictify, listify
from snappy_pipeline.workflows.abstract import BaseStep, BaseStepPart, LinkOutStepPart
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType
from snappy_pipeline.workflows.ngs_mapping import NgsMappingWorkflow, ResourceUsage
from snappy_pipeline.workflows.somatic_variant_annotation import (
    SomaticVariantAnnotationWorkflow,
)
from snappy_pipeline.workflows.somatic_variant_calling import (
    SomaticVariantCallingWorkflow,
)
from snappy_pipeline.workflows.somatic_variant_filtration import SomaticVariantFiltrationWorkflow

from .model import TumorMutationalBurden as TumorMutationalBurdenConfigModel

#: Extensions of files to create as main payload
EXT_VALUES = (".json", ".json.md5")

#: Names of the files to create for the extension
EXT_NAMES = ("json", "json_md5")

#: Default configuration for the tmb calculation step
DEFAULT_CONFIG = TumorMutationalBurdenConfigModel.default_config_yaml_string()


class TumorMutationalBurdenCalculationStepPart(BaseStepPart):
    """Calculation tumor mutational burden for each sample"""

    name = "tmb_gathering"

    actions = ("run",)

    def __init__(self, parent):
        super().__init__(parent)
        # Build shortcut from cancer bio sample name to matched cancer sample
        self.tumor_ngs_library_to_sample_pair = OrderedDict()
        for sheet in self.parent.shortcut_sheets:
            # update function of OrderedDict
            self.tumor_ngs_library_to_sample_pair.update(
                sheet.all_sample_pairs_by_tumor_dna_ngs_library
            )
        # Build mapping from donor name to donor.
        self.donors = OrderedDict()
        for sheet in self.parent.shortcut_sheets:
            for donor in sheet.donors:
                self.donors[donor.name] = donor

    @dictify
    def get_input_files(self, action):
        self._validate_action(action)

        base_name = "{tumor_library}"

        tpl = os.path.join("output", base_name, "out", base_name)

        key_ext = {"vcf": ".vcf.gz", "vcf_tbi": ".vcf.gz.tbi"}
        variant_path = self.parent.modules["somatic_variant"]
        for key, ext in key_ext.items():
            yield key, variant_path(tpl + ext)

    @dictify
    def get_output_files(self, action):
        # Validate action
        self._validate_action(action)

        base_name = "tmb.{tumor_library}"

        tpl = os.path.join("output", base_name, "out", base_name)

        key_ext = {"json": ".json"}
        for key, ext in key_ext.items():
            yield key, tpl + ext
            yield key + "_md5", tpl + ext + ".md5"

    @dictify
    def _get_log_file(self, action):
        self._validate_action(action)

        base_name = "tmb.{tumor_library}"

        tpl = os.path.join("output", base_name, "log", base_name)

        key_ext = (
            ("log", ".log"),
            ("conda_info", ".conda_info.txt"),
            ("conda_list", ".conda_list.txt"),
        )
        for key, ext in key_ext:
            yield key, tpl + ext

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
    def default_config_yaml(cls):
        """Return default config YAML, to be overwritten by project-specific one."""
        return DEFAULT_CONFIG

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
            previous_steps=(
                SomaticVariantCallingWorkflow,
                SomaticVariantAnnotationWorkflow,
                SomaticVariantFiltrationWorkflow,
                NgsMappingWorkflow,
            ),
            task_name=task_name,
            **kwargs,
        )
        # Register sub workflows
        config = self.config
        self.register_module("somatic_variant", str(config.somatic_variant_step))

        # Register sub step classes so the sub steps are available
        self.register_sub_step_classes((TumorMutationalBurdenCalculationStepPart, LinkOutStepPart))

    @listify
    def get_result_files(self):
        name_pattern = "tmb.{tumor_library.name}"

        yield from self._yield_result_files_matched(
            os.path.join("output", name_pattern, "out", name_pattern + "{ext}"),
            ext=EXT_VALUES,
        )
        yield from self._yield_result_files_matched(
            os.path.join("output", name_pattern, "log", name_pattern + "{ext}"),
            ext=(
                ".log",
                ".log.md5",
                ".conda_info.txt",
                ".conda_info.txt.md5",
                ".conda_list.txt",
                ".conda_list.txt.md5",
            ),
        )

    def _yield_result_files_matched(self, tpl, **kwargs):
        """Build output paths from path template and extension list.

        This function returns the results from the matched somatic variant callers such as
        Mutect.
        """
        for sheet in filter(is_not_background, self.shortcut_sheets):
            for sample_pair in sheet.all_sample_pairs:
                if (
                    not sample_pair.tumor_sample.dna_ngs_library
                    or not sample_pair.normal_sample.dna_ngs_library
                ):
                    msg = (
                        "INFO: sample pair for cancer bio sample {} has is missing primary"
                        "normal or primary cancer NGS library"
                    )
                    print(msg.format(sample_pair.tumor_sample.name), file=sys.stderr)
                    continue
                yield from expand(
                    tpl,
                    tumor_library=[sample_pair.tumor_sample.dna_ngs_library],
                    **kwargs,
                )
