from collections import OrderedDict
import os
import sys

from biomedsheets.shortcuts import CancerCaseSheet, CancerCaseSheetOptions, is_not_background
from snakemake.io import expand

from snappy_pipeline.utils import dictify, listify
from snappy_pipeline.workflows.abstract import BaseStep, BaseStepPart, LinkOutStepPart
from snappy_pipeline.workflows.ngs_mapping import NgsMappingWorkflow, ResourceUsage
from snappy_pipeline.workflows.somatic_variant_annotation import (
    ANNOTATION_TOOLS,
    SomaticVariantAnnotationWorkflow,
)
from snappy_pipeline.workflows.somatic_variant_calling import (
    SOMATIC_VARIANT_CALLERS_MATCHED,
    SomaticVariantCallingWorkflow,
)

#: Extensions of files to create as main payload
EXT_VALUES = (".json", ".json.md5")

#: Names of the files to create for the extension
EXT_NAMES = ("json", "json_md5")

#: Default configuration for the tmb calculation step
DEFAULT_CONFIG = r"""
step_config:
    tumor_mutational_burden:
        has_annotation: 'TRUE' # REQUIRED
        path_somatic_variant: ../somatic_variant_annotation   # REQUIRED
        tools_ngs_mapping: []      # default to those configured for ngs_mapping
        tools_somatic_variant_calling: []  # default to those configured for somatic_variant_calling
        tools_somatic_variant_annotation: [] # default to those configured for somatic_variant_annotation
        target_regions: # REQUIRED
        missense_regex: '.*[\|&]missense_variant[\|&].*' #change if the annotation tool doesn't use 'missense_variant' to indicate missense variant
"""


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
        # Adding part for runnng with annotation file instead of with variant calling file
        if self.w_config["step_config"]["tumor_mutational_burden"]["has_annotation"] == "TRUE":
            tpl = (
                "output/{mapper}.{var_caller}.{anno_tool}.{tumor_library}/out/"
                "{mapper}.{var_caller}.{anno_tool}.{tumor_library}"
            )
        else:
            tpl = (
                "output/{mapper}.{var_caller}.{tumor_library}/out/"
                "{mapper}.{var_caller}.{tumor_library}"
            )

        key_ext = {"vcf": ".vcf.gz", "vcf_tbi": ".vcf.gz.tbi"}
        # Adding part for runnng with annotation file instead of with variant calling file
        if self.w_config["step_config"]["tumor_mutational_burden"]["has_annotation"] == "TRUE":
            variant_path = self.parent.sub_workflows["somatic_variant_annotation"]
        else:
            variant_path = self.parent.sub_workflows["somatic_variant_calling"]
        for key, ext in key_ext.items():
            yield key, variant_path(tpl + ext)

    @dictify
    def get_output_files(self, action):
        # Validate action
        self._validate_action(action)
        if self.w_config["step_config"]["tumor_mutational_burden"]["has_annotation"] == "TRUE":
            prefix = (
                "work/{mapper}.{var_caller}.{anno_tool}.tmb.{tumor_library}/out/"
                "{mapper}.{var_caller}.{anno_tool}.tmb.{tumor_library}"
            )
        else:
            prefix = (
                "work/{mapper}.{var_caller}.tmb.{tumor_library}/out/"
                "{mapper}.{var_caller}.tmb.{tumor_library}"
            )
        key_ext = {"json": ".json"}
        for key, ext in key_ext.items():
            yield key, prefix + ext
            yield key + "_md5", prefix + ext + ".md5"

    @dictify
    def _get_log_file(self, action):
        self._validate_action(action)
        if self.w_config["step_config"]["tumor_mutational_burden"]["has_annotation"] == "TRUE":
            prefix = (
                "work/{mapper}.{var_caller}.{anno_tool}.tmb.{tumor_library}/log/"
                "{mapper}.{var_caller}.{anno_tool}.tmb.{tumor_library}"
            )
        else:
            prefix = (
                "work/{mapper}.{var_caller}.tmb.{tumor_library}/log/"
                "{mapper}.{var_caller}.tmb.{tumor_library}"
            )

        key_ext = (
            ("log", ".log"),
            ("conda_info", ".conda_info.txt"),
            ("conda_list", ".conda_list.txt"),
        )
        for key, ext in key_ext:
            yield key, prefix + ext

    def get_resource_usage(self, action):
        self._validate_action(action)
        mem_mb = 4 * 1024  # 4GB
        return ResourceUsage(
            threads=2,
            time="1:00:00",  # 1 hour
            memory=f"{mem_mb}M",
        )

    def get_params(self, action):
        self._validate_action(action)
        return getattr(self, "_get_params_run")

    def _get_params_run(self, wildcards):
        return {
            "missense_re": self.w_config["step_config"]["tumor_mutational_burden"]["missense_regex"]
        }


class TumorMutationalBurdenCalculationWorkflow(BaseStep):
    """Perform TMB calculation"""

    name = "tumormutation"
    sheet_shortcut_class = CancerCaseSheet
    sheet_shortcut_kwargs = {
        "options": CancerCaseSheetOptions(allow_missing_normal=True, allow_missing_tumor=True)
    }

    @classmethod
    def default_config_yaml(cls):
        """Return default config YAML, to be overwritten by project-specific one."""
        return DEFAULT_CONFIG

    def __init__(self, workflow, config, config_lookup_paths, config_paths, workdir):
        super().__init__(
            workflow,
            config,
            config_lookup_paths,
            config_paths,
            workdir,
            (SomaticVariantCallingWorkflow, SomaticVariantAnnotationWorkflow, NgsMappingWorkflow),
        )
        # Register sub step classes so the sub steps are available
        self.register_sub_step_classes((TumorMutationalBurdenCalculationStepPart, LinkOutStepPart))
        # Register sub workflows
        if self.w_config["step_config"]["tumor_mutational_burden"]["has_annotation"] == "TRUE":
            self.register_sub_workflow(
                "somatic_variant_annotation",
                self.w_config["step_config"]["tumor_mutational_burden"]["path_somatic_variant"],
            )
            if not self.w_config["step_config"]["tumor_mutational_burden"][
                "tools_somatic_variant_annotation"
            ]:
                self.w_config["step_config"]["tumor_mutational_burden"][
                    "tools_somatic_variant_annotation"
                ] = self.w_config["step_config"]["somatic_variant_annotation"]["tools"]
        else:
            self.register_sub_workflow(
                "somatic_variant_calling",
                self.w_config["step_config"]["tumor_mutational_burden"]["path_somatic_variant"],
            )
        # Copy over "tools" setting from somatic_variant_calling/ngs_mapping if not set here
        if not self.w_config["step_config"]["tumor_mutational_burden"]["tools_ngs_mapping"]:
            self.w_config["step_config"]["tumor_mutational_burden"][
                "tools_ngs_mapping"
            ] = self.w_config["step_config"]["ngs_mapping"]["tools"]["dna"]
        if not self.w_config["step_config"]["tumor_mutational_burden"][
            "tools_somatic_variant_calling"
        ]:
            self.w_config["step_config"]["tumor_mutational_burden"][
                "tools_somatic_variant_calling"
            ] = self.w_config["step_config"]["somatic_variant_calling"]["tools"]

    @listify
    def get_result_files(self):
        callers = set(
            self.w_config["step_config"]["tumor_mutational_burden"]["tools_somatic_variant_calling"]
        )
        if self.w_config["step_config"]["tumor_mutational_burden"]["has_annotation"] == "TRUE":
            anno_callers = set(
                self.w_config["step_config"]["tumor_mutational_burden"][
                    "tools_somatic_variant_annotation"
                ]
            )
            name_pattern = "{mapper}.{caller}.{anno_caller}.tmb.{tumor_library.name}"
            yield from self._yield_result_files_matched(
                os.path.join("output", name_pattern, "out", name_pattern + "{ext}"),
                mapper=self.w_config["step_config"]["tumor_mutational_burden"]["tools_ngs_mapping"],
                caller=callers & set(SOMATIC_VARIANT_CALLERS_MATCHED),
                anno_caller=anno_callers & set(ANNOTATION_TOOLS),
                ext=EXT_VALUES,
            )
            yield from self._yield_result_files_matched(
                os.path.join("output", name_pattern, "log", name_pattern + "{ext}"),
                mapper=self.w_config["step_config"]["tumor_mutational_burden"]["tools_ngs_mapping"],
                caller=callers & set(SOMATIC_VARIANT_CALLERS_MATCHED),
                anno_caller=anno_callers & set(ANNOTATION_TOOLS),
                ext=(
                    ".log",
                    ".log.md5",
                    ".conda_info.txt",
                    ".conda_info.txt.md5",
                    ".conda_list.txt",
                    ".conda_list.txt.md5",
                ),
            )
        else:
            name_pattern = "{mapper}.{caller}.tmb.{tumor_library.name}"
            yield from self._yield_result_files_matched(
                os.path.join("output", name_pattern, "out", name_pattern + "{ext}"),
                mapper=self.w_config["step_config"]["tumor_mutational_burden"]["tools_ngs_mapping"],
                caller=callers & set(SOMATIC_VARIANT_CALLERS_MATCHED),
                ext=EXT_VALUES,
            )
            yield from self._yield_result_files_matched(
                os.path.join("output", name_pattern, "log", name_pattern + "{ext}"),
                mapper=self.w_config["step_config"]["tumor_mutational_burden"]["tools_ngs_mapping"],
                caller=callers & set(SOMATIC_VARIANT_CALLERS_MATCHED),
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

    def check_config(self):
        """Check that the path to the NGS mapping is present"""
        self.ensure_w_config(
            ("step_config", "tumor_mutational_burden", "path_somatic_variant"),
            "Path to variant (directory of vcf files) not configured but required for tmb calculation",
        )

        self.ensure_w_config(
            ("step_config", "tumor_mutational_burden", "target_regions"),
            "Path to target_regions file (bed format)"
            "not configured but required for tmb calculation",
        )

        self.ensure_w_config(
            ("step_config", "tumor_mutational_burden", "has_annotation"),
            "TMB needs to know wether the vcf is annotated or not",
        )
